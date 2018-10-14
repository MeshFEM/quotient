/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_MINIMUM_DEGREE_IMPL_H_
#define QUOTIENT_MINIMUM_DEGREE_IMPL_H_

#include <iostream>
#include <unordered_map>
#include <vector>

#include "quotient/config.hpp"
#include "quotient/coordinate_graph.hpp"
#include "quotient/external_degree.hpp"
#include "quotient/minimum_degree.hpp"
#include "quotient/quotient_graph.hpp"
#include "quotient/timer.hpp"

namespace quotient {

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

inline double
MinimumDegreeAnalysis::FractionOfPivotsWithMultipleElements() const {
  Int num_pivots_with_multiple_elements = 0;
  for (const Int& pivot_element_list_size : pivot_element_list_sizes) {
    if (pivot_element_list_size > 2) {
      ++num_pivots_with_multiple_elements;
    }
  }
  return num_pivots_with_multiple_elements /
      (1. * pivot_element_list_sizes.size());
}

inline double
MinimumDegreeAnalysis::FractionOfDegreeUpdatesWithMultipleElements() const {
  const Int num_total_degree_updates =
      num_degree_updates_with_multiple_elements +
      num_degree_updates_without_multiple_elements;
  return num_degree_updates_with_multiple_elements /
      (1. * num_total_degree_updates);
}

namespace minimum_degree {

namespace compute_pivot_structure {

// Computes the structure of the pivot:
//
//   L_p := (A_p \cup (\cup_{e in E_p} L_e)) \ supernode(p).
//
// in the case where there are modest numbers of members of E_p.
inline void FewElements(Int pivot, QuotientGraph* graph) {
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
inline void ManyElements(Int pivot, QuotientGraph* graph) {
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
        graph->structures[pivot].begin() + offsets[index]);
  }

  // Sort the non-unique list of elements and then erase duplicates.
  //
  // TODO(Jack Poulson): Use a multithreaded sort.
  std::sort(graph->structures[pivot].begin(), graph->structures[pivot].end());
  EraseDuplicatesInSortedVector(&graph->structures[pivot]);
}

}  // namespace compute_pivot_structure

// Computes the structure of the pivot:
//
//   L_p := (A_p \cup (\cup_{e in E_p} L_e)) \ supernode(p).
//
inline void ComputePivotStructure(Int pivot, QuotientGraph* graph) {
  const Int kElementThreshold = 4;
  if (graph->element_lists[pivot].size() < kElementThreshold) {
    compute_pivot_structure::FewElements(pivot, graph);
    return;
  }
  compute_pivot_structure::ManyElements(pivot, graph);
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
    const std::vector<Int>& aggressive_absorption_elements,
    bool store_aggressive_absorptions,
    QuotientGraph* graph) {
  for (const Int& i : supernodal_pivot_structure) {
    // Element absorption:
    //   E_i := (E_i \ E_p) \cup {p}.
    FilterSetInPlace(graph->element_lists[pivot], &graph->element_lists[i]);
    InsertNewEntryIntoSet(pivot, &graph->element_lists[i]);
  }

  for (const Int& element : aggressive_absorption_elements) {
    // Aggressive element absorption:
    //   E_p := (E_p \ E_e) \cup {e}.
    if (store_aggressive_absorptions) {
      graph->aggressive_absorptions.emplace_back(pivot, element);
    }
    FilterSetInPlace(
        graph->element_lists[element], &graph->element_lists[pivot]);
    InsertEntryIntoSet(element, &graph->element_lists[pivot]);
  }
}

// Compute the external degree approximations of the supernodes
// adjacent to the current pivot.
inline void UpdateExternalDegreeApproximations(
    Int pivot,
    const std::vector<Int>& supernodal_pivot_structure,
    const std::vector<Int>& external_structure_sizes,
    ExternalDegreeType degree_type,
    QuotientGraph* graph) {
  // Compute and store the external degrees.
  // TODO(Jack Poulson): Parallelize this loop with OpenMP.
  const std::size_t struct_size = supernodal_pivot_structure.size();
  std::vector<Int> external_degrees(struct_size);
  for (std::size_t index = 0; index < struct_size; ++index) {
    const Int i = supernodal_pivot_structure[index];

    // Compute the external degree (or approximation) of supervariable i:
    //   d_i := |A_i \ supernode(i)| +
    //       |(\cup_{e in E_i} L_e) \ supernode(i)|.
    external_degrees[index] = ExternalDegree(
        *graph, i, pivot, external_structure_sizes, degree_type);
  }

  // Insert the external degrees into the heap.
  //
  // TODO(Jack Poulson): Add support for a batch interface to
  // 'external_degree_heap.SetValue' so that comparison propagations can be
  // potentially shared.
  for (std::size_t index = 0; index < struct_size; ++index) {
    const Int i = supernodal_pivot_structure[index];
    graph->external_degree_heap.SetValue(i, external_degrees[index]);
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
  const std::size_t struct_size = supernodal_pivot_structure.size();
  std::vector<Int> bucket_list(struct_size);
  std::unordered_map<Int, std::vector<Int>> variable_hash_map;
  for (std::size_t i_index = 0; i_index < struct_size; ++i_index) {
    const Int i = supernodal_pivot_structure[i_index];

    // Append this principal variable to its hash bucket.
    // TODO(Jack Poulson): Add configurable support for other hashes.
    // For example, one could mod by (std::numeric_limits<Int>::max() - 1).
    const std::size_t bucket_key = graph->AshcraftVariableHash(i);
    bucket_list[i_index] = bucket_key;
    variable_hash_map[bucket_key].push_back(i_index);
  }

  // TODO(Jack Poulson): Refactor and parallelize this loop with OpenMP.
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

}  // namespace minimum_degree

inline MinimumDegreeAnalysis MinimumDegree(
  const CoordinateGraph& graph, const MinimumDegreeControl& control) {
  using namespace minimum_degree;
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
  constexpr char kDetectAndMergeVariables[] = "DetectAndMergeVariables";
  if (control.time_stages) {
    timers[kComputePivotStructure].Reset();
    timers[kUpdateAdjacencyLists].Reset();
    timers[kExternalStructureSizes].Reset();
    timers[kUpdateElementLists].Reset();
    timers[kUpdateExternalDegrees].Reset();
    timers[kDetectAndMergeVariables].Reset();
  }
  const bool compute_external_structure_sizes =
      control.aggressive_absorption ||
      control.degree_type == kAmestoyExternalDegree;
  std::vector<Int> external_structure_sizes;
  if (compute_external_structure_sizes) {
    if (control.time_stages) timers[kExternalStructureSizes].Start();
    quotient_graph.InitializeExternalStructureSizes(&external_structure_sizes);
    if (control.time_stages) timers[kExternalStructureSizes].Stop();
  }
  std::vector<Int> aggressive_absorption_elements;
  if (control.store_pivot_element_list_sizes) {
    analysis.pivot_element_list_sizes.reserve(num_orig_vertices);
  }
  if (control.store_num_degree_updates_with_multiple_elements) {
    analysis.num_degree_updates_with_multiple_elements = 0;
    analysis.num_degree_updates_without_multiple_elements = 0;
  }
  while (quotient_graph.num_eliminated_vertices < num_orig_vertices) {
    // Retrieve a variable with minimal (approximate) external degree.
    const std::pair<Int, Int> pivot_pair =
        quotient_graph.external_degree_heap.MinimalEntry(); 
    const Int pivot = pivot_pair.first;
    if (control.store_pivot_element_list_sizes) {
      analysis.pivot_element_list_sizes.push_back(
          quotient_graph.element_lists[pivot].size());
    }

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
    if (compute_external_structure_sizes) {
      if (control.time_stages) timers[kExternalStructureSizes].Start();
      quotient_graph.ExternalStructureSizes(
          supernodal_pivot_structure, control.aggressive_absorption,
          &external_structure_sizes, &aggressive_absorption_elements);
      if (control.time_stages) timers[kExternalStructureSizes].Stop();
    }

    if (control.time_stages) timers[kUpdateElementLists].Start();
    UpdateElementListsAfterSelectingPivot(
        pivot, supernodal_pivot_structure, aggressive_absorption_elements,
        control.store_aggressive_absorptions, &quotient_graph);
    if (control.time_stages) timers[kUpdateElementLists].Stop();

    if (control.time_stages) timers[kUpdateExternalDegrees].Start();
    UpdateExternalDegreeApproximations(
        pivot, supernodal_pivot_structure, external_structure_sizes,
        control.degree_type, &quotient_graph);
    if (control.time_stages) timers[kUpdateExternalDegrees].Stop();
    if (control.store_num_degree_updates_with_multiple_elements) {
      for (const Int& i : supernodal_pivot_structure) {
        if (quotient_graph.element_lists[i].size() > 2) {
          ++analysis.num_degree_updates_with_multiple_elements;
        } else {
          ++analysis.num_degree_updates_without_multiple_elements;
        }
      }
    }

    if (control.allow_supernodes) {
      if (control.time_stages) timers[kDetectAndMergeVariables].Start();
      DetectAndMergeVariables(
          supernodal_pivot_structure, control.store_variable_merges,
          &quotient_graph);
      if (control.time_stages) timers[kDetectAndMergeVariables].Stop();
    }

    ConvertPivotIntoElement(pivot, &quotient_graph);

    analysis.elimination_order.push_back(pivot);

    if (compute_external_structure_sizes) {
      if (control.time_stages) timers[kExternalStructureSizes].Start();
      quotient_graph.ResetExternalStructureSizes(
          supernodal_pivot_structure, &external_structure_sizes);
      if (control.time_stages) timers[kExternalStructureSizes].Stop();
    }
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

#endif // ifndef QUOTIENT_MINIMUM_DEGREE_IMPL_H_
