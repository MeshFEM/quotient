/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_MINIMUM_DEGREE_H_
#define QUOTIENT_MINIMUM_DEGREE_H_

#include <iostream>
#include <unordered_map>
#include <vector>

#include "quotient/config.hpp"
#include "quotient/coordinate_graph.hpp"
#include "quotient/external_degree.hpp"
#include "quotient/quotient_graph.hpp"
#include "quotient/random_access_heap.hpp"
#include "quotient/timer.hpp"

namespace quotient {

// A data structure for controlling the MinimumDegree reordering routine.
struct MinimumDegreeControl {
  // The type of approximation to use for the external degree estimates.
  ExternalDegreeType degree_type = kAmestoyExternalDegree;

  // Whether nontrivial supernodes are allowed. It is highly recommended that
  // this remain true.
  bool allow_supernodes = true;

  // Whether aggressive element absorptions are allowed.
  bool aggressive_absorption = true;

  // Whether the list of pairs of aggressive element absorptions should be
  // returned in the MinimumDegreeAnalysis result of MinimumDegree.
  bool store_aggressive_absorptions = false;

  // Whether the list of pairs of variable merges should be returned in the
  // MinimumDegreeAnalysis result of MinimumDegree.
  bool store_variable_merges = false;

  // Whether a breakdown of the elapsed seconds of each stage of the reordering
  // should be saved.
  bool time_stages = false;
};


// The result of running the MinimumDegree reordering algorithm. It contains
// the ordered list of eliminated principal vertices, the list of supernodes,
// and the supernodal nonzero structure of each principal column.
struct MinimumDegreeAnalysis {
  // The recommended elimination order of the supernodes (with each supernode
  // represented by its principal member).
  std::vector<Int> elimination_order;

  // The list of supernodes. Entry 'i' is nonempty if and only if 'i' is the
  // principal member of the supernode, in which case 'supernodes[i]' contains
  // the sorted list of members of the supernode.
  std::vector<std::vector<Int>> supernodes;

  // The structures of the supernodes. Entry 'index' corresponds to the
  // structure of supernode 'elimination_order[index]'.
  std::vector<std::vector<Int>> principal_structures;

  // An optional list (based on the value of
  // 'MinimumDegreeControl.store_aggressive_absorptions') of aggressive element
  // absorption pairs: each pair (e, f) consists of the absorbing element, e,
  // and the absorbed element, f.
  std::vector<std::pair<Int, Int>> aggressive_absorptions;

  // An optional list (based on the value of
  // 'MinimumDegreeControl.store_variable_merges') of supervariable merge
  // pairs: each pair (i, j) consists of the absorbing supervariable, i, and
  // the absorbed supervariable, j.
  std::vector<std::pair<Int, Int>> variable_merges;

  // We will push to elimination order as the reordering algorithm progresses,
  // so we will allocate an upper bound for the amount of required space.
  // The 'supernodes' and 'structures' variables will be copied over from
  // the quotient graph just before the analysis completes.
  MinimumDegreeAnalysis(Int num_vertices);

  // Returns the number of structural nonzeros in the strictly lower-triangular
  // factor.
  Int NumStrictlyLowerNonzeros() const;

  // Returns the principal member of the largest supernode.
  Int LargestSupernode() const;

  // Returns the number of members of the largest supernode.
  Int LargestSupernodeSize() const;

  // An optional (based upon the value of 'MinimumDegreeControl.time_stages')
  // map from the stage names to the corresponding elapsed seconds.
  std::unordered_map<std::string, double> elapsed_seconds;
};

inline MinimumDegreeAnalysis::MinimumDegreeAnalysis(Int num_vertices) {
  elimination_order.reserve(num_vertices);
}

inline Int MinimumDegreeAnalysis::NumStrictlyLowerNonzeros() const {
  Int num_nonzeros = 0;
  for (std::size_t index = 0; index < elimination_order.size(); ++index) {
    const Int j = elimination_order[index];
    const Int supernode_j_size = supernodes[j].size();

    // Add the strictly-lower triangular portion of the diagonal block.
    num_nonzeros += ((supernode_j_size - 1) * supernode_j_size) / 2;

    // Add the rectangular portion below the diagonal block.
    num_nonzeros += supernode_j_size * principal_structures[index].size();
  }
  return num_nonzeros;
}

inline Int MinimumDegreeAnalysis::LargestSupernode() const {
  Int largest_supernode = -1;
  std::size_t largest_supernode_size = 0;
  for (std::size_t i = 0; i < supernodes.size(); ++i) {
    if (supernodes[i].size() > largest_supernode_size) {
      largest_supernode = i;
      largest_supernode_size = supernodes[i].size();
    }
  }
  return largest_supernode;
}

inline Int MinimumDegreeAnalysis::LargestSupernodeSize() const {
  return supernodes[LargestSupernode()].size();
}

// Computes the structure of the pivot:
//
//   L_p := (A_p \cup (\cup_{e in E_p} L_e)) \ supernode(p).
//
// in the case where there are modest numbers of members of E_p.
inline void ComputePivotStructureFewElements(Int pivot, QuotientGraph* graph) {
  FilterSet(
      graph->adjacency_lists[pivot],
      graph->supernodes[pivot],
      &graph->structures[pivot]);

  std::vector<Int> temp0, temp1; 
  for (const Int& element : graph->element_lists[pivot]) {
    FilterSet(graph->structures[element], graph->supernodes[pivot], &temp0);
    temp1 = graph->structures[pivot];
    MergeSets(temp1, temp0, &graph->structures[pivot]);
  }
}

// Computes the structure of the pivot:
//
//   L_p := (A_p \cup (\cup_{e in E_p} L_e)) \ supernode(p).
//
// in the case where there are many numbers of members of E_p.
inline void ComputePivotStructureManyElements(Int pivot, QuotientGraph* graph) {
  // Set up the list of index sets we will be merging (after filtering
  // supernode(pivot)).
  std::vector<std::vector<Int>*> index_sets;
  index_sets.reserve(graph->element_lists[pivot].size() + 1);
  index_sets.push_back(&graph->adjacency_lists[pivot]);
  for (const Int& element : graph->element_lists[pivot]) {
    index_sets.push_back(&graph->structures[element]);
  }

  // Count the number of (non-unique) filtered entries and compute the scan
  // of the offsets.
  Int num_non_unique = 0;
  std::vector<Int> offsets;
  offsets.reserve(index_sets.size());
  for (const std::vector<Int>* index_set : index_sets) {
    offsets.push_back(num_non_unique);
    num_non_unique += SizeOfDifference(*index_set, graph->supernodes[pivot]);
  }

  // Fill the unsorted, non-unique list of elements.
  //
  // TODO(Jack Poulson): Use OpenMP to parallelize this loop.
  graph->structures[pivot].clear();
  graph->structures[pivot].resize(num_non_unique);
  for (std::size_t index = 0; index < index_sets.size(); ++index) {
    std::set_difference(
        index_sets[index]->begin(), index_sets[index]->end(),
        graph->supernodes[pivot].begin(), graph->supernodes[pivot].end(),
        graph->structures[pivot].data() + offsets[index]);
  }

  // Sort the non-unique list of elements and then erase duplicates.
  //
  // TODO(Jack Poulson): Use a multithreaded sort.
  std::sort(graph->structures[pivot].begin(), graph->structures[pivot].end());
  EraseDuplicatesInSortedVector(&graph->structures[pivot]);
}

// Computes the structure of the pivot:
//
//   L_p := (A_p \cup (\cup_{e in E_p} L_e)) \ supernode(p).
//
inline void ComputePivotStructure(Int pivot, QuotientGraph* graph) {
  const Int kElementThreshold = 4;
  if (graph->element_lists[pivot].size() < kElementThreshold) {
    ComputePivotStructureFewElements(pivot, graph);
    return;
  }
  ComputePivotStructureManyElements(pivot, graph);
}

// Update the adjacency lists after computing the supernodal pivot structure.
inline void UpdateAdjacencyListsAfterSelectingPivot(
    Int pivot,
    const std::vector<Int>& supernodal_pivot_structure,
    QuotientGraph* graph) {
  for (const Int& i : supernodal_pivot_structure) {
    // Remove redundant adjacency entries:
    //   A_i := (A_i \ L_p) \ supernode(p).
    DoubleFilterSetInPlace(
      graph->structures[pivot],
      graph->supernodes[pivot],
      &graph->adjacency_lists[i]);
  }
}

// Update the element lists after computing the supernodal pivot structure.
inline void UpdateElementListsAfterSelectingPivot(
    Int pivot,
    const std::vector<Int>& supernodal_pivot_structure,
    const std::unordered_map<Int, Int>& external_structure_sizes,
    bool aggressive_absorption,
    bool store_aggressive_absorptions,
    QuotientGraph* graph) {
  for (const Int& i : supernodal_pivot_structure) {
    // Element absorption:
    //   E_i := (E_i \ E_p) \cup {p}.
    FilterSetInPlace(graph->element_lists[pivot], &graph->element_lists[i]);
    InsertNewEntryIntoSet(pivot, &graph->element_lists[i]);
  }

  if (aggressive_absorption) {
    // Follow the advice at the beginning of Section 5 of [AMD-96] and absorb
    // any element e that satisfies |L_e \ L_p| = 0.
    for (const std::pair<Int, Int>& pairing : external_structure_sizes) {
      if (pairing.second != 0) {
        continue;
      }

      // Aggressive element absorption:
      //   E_p := (E_p \ E_e) \cup {e}.
      const Int element = pairing.first;
      if (store_aggressive_absorptions) {
        graph->aggressive_absorptions.emplace_back(pivot, element);
      }
      FilterSetInPlace(
          graph->element_lists[element], &graph->element_lists[pivot]);
      InsertEntryIntoSet(element, &graph->element_lists[pivot]);
    }
  }
}

// Compute the external degree approximations of the supernodes
// adjacent to the current pivot.
//
// TODO(Jack Poulson): Add support for a batch interface to
// 'external_degree_heap.SetValue' so that comparison propagations can be
// potentially shared.
inline void UpdateExternalDegreeApproximations(
    Int pivot,
    const std::vector<Int>& supernodal_pivot_structure,
    const std::unordered_map<Int, Int>& external_structure_sizes,
    ExternalDegreeType degree_type,
    QuotientGraph* graph) {
  for (const Int& i : supernodal_pivot_structure) {
    // Compute the external degree (or approximation) of supervariable i:
    //   d_i := |A_i \ supernode(i)| +
    //       |(\cup_{e in E_i} L_e) \ supernode(i)|.
    const Int external_degree = ExternalDegree(
        *graph, i, pivot, external_structure_sizes, degree_type);
    graph->external_degree_heap.SetValue(i, external_degree);
  }
}

// Detects and merges pairs of supervariables in the pivot structure who are
// indistinguishable with respect to the quotient graph.
//
// While the supernodal merges will potentially shrink the supernodal adjacency
// lists (and thus change the associated Ashcraft hash of variables), if two
// variables are indistinguishable, their cached bucket might be wrong, but
// they would be wrong together.
//
// The test for indistinguishability does not depend upon the variable
// supernodal structure and is thus invariant to supervariable merges.
inline void DetectAndMergeVariables(
    const std::vector<Int>& supernodal_pivot_structure,
    bool store_variable_merges,
    QuotientGraph* graph) {
  // Fill a set of buckets for the hashes of the supernodes adjacent to
  // the current pivot.
  const Int supernodal_pivot_struct_size = supernodal_pivot_structure.size();
  std::vector<Int> bucket_list(supernodal_pivot_struct_size);
  std::unordered_map<Int, std::vector<Int>> variable_hash_map;
  for (Int i_index = 0; i_index < supernodal_pivot_struct_size; ++i_index) {
    const Int i = supernodal_pivot_structure[i_index];

    // Append this principal variable to its hash bucket.
    // TODO(Jack Poulson): Add configurable support for other hashes.
    // For example, one could mod by (std::numeric_limits<Int>::max() - 1).
    const std::size_t bucket_key = graph->AshcraftVariableHash(i);
    bucket_list[i_index] = bucket_key;
    variable_hash_map[bucket_key].push_back(i_index);
  }

  std::vector<Int> temp;
  for (const std::pair<Int, std::vector<Int>>& entry : variable_hash_map) {
    const std::vector<Int>& bucket = entry.second;
    std::vector<bool> merged_supernode(bucket.size(), false);
    for (std::size_t i_bucket  = 0; i_bucket < bucket.size(); ++i_bucket) {
      if (merged_supernode[i_bucket]) {
        continue;
      }
      const Int i_index = bucket[i_bucket];
      const Int i = supernodal_pivot_structure[i_index];
      for (std::size_t j_bucket = 0; j_bucket < bucket.size(); ++j_bucket) {
        const Int j_index = bucket[j_bucket];
        // Avoid processing the same pair twice, and skip any already-merged
        // supernodes.
        if (j_index <= i_index || merged_supernode[j_bucket]) {
          continue;
        }
        const Int j = supernodal_pivot_structure[j_index];
        if (graph->StructuralSupervariablesAreQuotientIndistinguishable(i, j)) {
          // Absorb supernode(j) into supernode(i). 
          if (store_variable_merges) {
            graph->variable_merges.emplace_back(i, j);
          }
          temp = graph->supernodes[i];
          MergeSets(temp, graph->supernodes[j], &graph->supernodes[i]);
          graph->external_degree_heap.UpdateValue(
              i, -graph->supernodes[j].size());
          graph->external_degree_heap.DisableIndex(j);
          SwapClearVector(&graph->supernodes[j]);
          SwapClearVector(&graph->adjacency_lists[j]);
          SwapClearVector(&graph->element_lists[j]);
          merged_supernode[j_bucket] = true;
        }
      }
    }
  }
}

// Converts the 'pivot' (super)variable into an element.
inline void ConvertPivotIntoElement(Int pivot, QuotientGraph* graph) {
  graph->external_degree_heap.DisableIndex(pivot);
  SwapClearVector(&graph->adjacency_lists[pivot]);
  SwapClearVector(&graph->element_lists[pivot]);
  graph->num_eliminated_vertices += graph->supernodes[pivot].size();
}

// Returns a supernodal reordering and the corresponding supernodal nonzero
// structure of the implied factorization using the (Approximate) Minimum
// Degree reordering algorithm. Please see [ADD-96] for details.
//
// The input graph must be explicitly symmetric.
//
MinimumDegreeAnalysis MinimumDegree(
  const CoordinateGraph& graph, const MinimumDegreeControl& control) {
#ifdef QUOTIENT_DEBUG
  if (graph.NumSources() != graph.NumTargets()) {
    std::cerr << "ERROR: MinimumDegree requires a symmetric input graph."
              << std::endl;
    return MinimumDegreeAnalysis(0);
  }
#endif
  QuotientGraph quotient_graph(graph);
  const Int num_orig_vertices = quotient_graph.num_original_vertices;

  // Eliminate the variables.
  MinimumDegreeAnalysis analysis(num_orig_vertices);
  std::unordered_map<std::string, Timer> timers;
  constexpr char kComputePivotStructure[] = "ComputePivotStructure";
  constexpr char kUpdateAdjacencyLists[] = "UpdateAdjacencyLists";
  constexpr char kExternalStructureSizes[] = "ExternalStructureSizes";
  constexpr char kUpdateElementLists[] = "UpdateElementLists";
  constexpr char kUpdateExternalDegrees[] = "UpdateExternalDegrees";
  constexpr char kDetectAndMergeVariables[] = "UpdateExternalVariables";
  if (control.time_stages) {
    timers[kComputePivotStructure].Reset();
    timers[kUpdateAdjacencyLists].Reset();
    timers[kExternalStructureSizes].Reset();
    timers[kUpdateElementLists].Reset();
    timers[kUpdateExternalDegrees].Reset();
    timers[kDetectAndMergeVariables].Reset();
  }
  while (quotient_graph.num_eliminated_vertices < num_orig_vertices) {
    // Retrieve a variable with minimal (approximate) external degree.
    const std::pair<Int, Int> pivot_pair =
        quotient_graph.external_degree_heap.MinimalEntry(); 
    const Int pivot = pivot_pair.first;

    if (control.time_stages) timers[kComputePivotStructure].Start();
    ComputePivotStructure(pivot, &quotient_graph);
    if (control.time_stages) timers[kComputePivotStructure].Stop();
    const std::vector<Int> supernodal_pivot_structure =
        quotient_graph.FormSupernodalStructure(pivot);

    if (control.time_stages) timers[kUpdateAdjacencyLists].Start();
    UpdateAdjacencyListsAfterSelectingPivot(
        pivot, supernodal_pivot_structure, &quotient_graph);
    if (control.time_stages) timers[kUpdateAdjacencyLists].Stop();

    // Compute the external structure cardinalities of the elements.
    // (but only if the Amestoy external degree approximation is requested).
    std::unordered_map<Int, Int> external_structure_sizes;
    if (control.aggressive_absorption ||
        control.degree_type == kAmestoyExternalDegree) {
      if (control.time_stages) timers[kExternalStructureSizes].Start();
      external_structure_sizes = quotient_graph.ExternalStructureSizes(
          supernodal_pivot_structure);
      if (control.time_stages) timers[kExternalStructureSizes].Stop();
    }

    if (control.time_stages) timers[kUpdateElementLists].Start();
    UpdateElementListsAfterSelectingPivot(
        pivot, supernodal_pivot_structure, external_structure_sizes,
        control.aggressive_absorption, control.store_aggressive_absorptions,
        &quotient_graph);
    if (control.time_stages) timers[kUpdateElementLists].Stop();

    if (control.time_stages) timers[kUpdateExternalDegrees].Start();
    UpdateExternalDegreeApproximations(
        pivot, supernodal_pivot_structure, external_structure_sizes,
        control.degree_type, &quotient_graph);
    if (control.time_stages) timers[kUpdateExternalDegrees].Stop();

    if (control.allow_supernodes) {
      if (control.time_stages) timers[kDetectAndMergeVariables].Start();
      DetectAndMergeVariables(
          supernodal_pivot_structure, control.store_variable_merges,
          &quotient_graph);
      if (control.time_stages) timers[kDetectAndMergeVariables].Stop();
    }

    ConvertPivotIntoElement(pivot, &quotient_graph);

    analysis.elimination_order.push_back(pivot);
  }

  // Extract the relevant information from the QuotientGraph.
  analysis.supernodes = quotient_graph.supernodes;
  analysis.principal_structures.resize(analysis.elimination_order.size());
  for (std::size_t index = 0; index < analysis.elimination_order.size();
       ++index) {
    analysis.principal_structures[index] =
        quotient_graph.structures[analysis.elimination_order[index]];
  }
  analysis.aggressive_absorptions = quotient_graph.aggressive_absorptions;
  analysis.variable_merges = quotient_graph.variable_merges;
  for (const std::pair<std::string, Timer>& pairing : timers) {
    analysis.elapsed_seconds[pairing.first] = pairing.second.TotalSeconds();
  }

  return analysis;
}

} // namespace quotient

#endif // ifndef QUOTIENT_MINIMUM_DEGREE_H_
