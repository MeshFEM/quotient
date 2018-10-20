/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_MINIMUM_DEGREE_IMPL_H_
#define QUOTIENT_MINIMUM_DEGREE_IMPL_H_

#ifdef _OPENMP
#include "omp.h"
#endif

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
  graph->structures[pivot].clear();
  graph->structures[pivot].resize(num_non_unique);
  #pragma omp parallel for schedule(dynamic)
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
  const std::size_t struct_size = supernodal_pivot_structure.size();
  std::vector<Int> external_degrees(struct_size);
  #pragma omp parallel for schedule(dynamic)
  for (std::size_t index = 0; index < struct_size; ++index) {
    const Int i = supernodal_pivot_structure[index];

    // Compute the external degree (or approximation) of supervariable i:
    //   d_i := |A_i \ supernode(i)| +
    //       |(\cup_{e in E_i} L_e) \ supernode(i)|.
    external_degrees[index] = ExternalDegree(
        *graph, i, pivot, external_structure_sizes, degree_type);
  }

  // Insert the external degrees into the heap.
  graph->external_degree_heap.ReserveValueChanges(struct_size);
  for (std::size_t index = 0; index < struct_size; ++index) {
    graph->external_degree_heap.QueueValueAssignment(
        supernodal_pivot_structure[index], external_degrees[index]);
  }
  graph->external_degree_heap.FlushValueChangeQueue();
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
    VariableHashType hash_type,
    QuotientGraph* graph) {
  // Compute the hashes for each variable.
  const std::size_t struct_size = supernodal_pivot_structure.size();
  std::vector<std::size_t> bucket_keys(struct_size);
  #pragma omp parallel for schedule(dynamic)
  for (std::size_t i_index = 0; i_index < struct_size; ++i_index) {
    const Int i = supernodal_pivot_structure[i_index];
    bucket_keys[i_index] = graph->VariableHash(i, hash_type);
  }

  // Fill a set of buckets for the hashes of the supernodes adjacent to
  // the current pivot.
  std::vector<std::vector<Int>> buckets;
  buckets.reserve(struct_size);
  std::unordered_map<std::size_t, std::size_t> key_to_bucket_index;
  for (std::size_t i_index = 0; i_index < struct_size; ++i_index) {
    const std::size_t bucket_key = bucket_keys[i_index];
    auto iter = key_to_bucket_index.find(bucket_key);
    if (iter == key_to_bucket_index.end()) {
      key_to_bucket_index[bucket_key] = buckets.size();
      buckets.emplace_back(1, i_index);
    } else {
      buckets[iter->second].push_back(i_index);
    }
  }

  // A data structure for representing the relevant information of a merge of
  // one supervariable (the 'absorbed_index') into another
  // (the 'primary_index').
  struct VariableMergeInfo {
    // The principal index of the supervariable that is increasing in size.
    Int primary_index;

    // The principal index of the supervariable being absorbed.
    Int absorbed_index;

    // The number of members of the absorbed supervarible.
    Int absorbed_size;

    VariableMergeInfo(
        Int primary_index_value,
        Int absorbed_index_value,
        Int absorbed_size_value) :
    primary_index(primary_index_value),
    absorbed_index(absorbed_index_value),
    absorbed_size(absorbed_size_value) {}
  };

  std::vector<Int> temp;
  std::vector<VariableMergeInfo> variable_merges;
  #pragma omp parallel for schedule(dynamic) private(temp)
  for (std::size_t index = 0; index < buckets.size(); ++index) {
    const std::vector<Int>& bucket = buckets[index];
    const Int bucket_size = bucket.size();
    std::vector<bool> merged_supernode(bucket_size, false);
    for (Int i_bucket  = 0; i_bucket < bucket_size; ++i_bucket) {
      if (merged_supernode[i_bucket]) {
        continue;
      }
      const Int i_index = bucket[i_bucket];
      const Int i = supernodal_pivot_structure[i_index];
      for (Int j_bucket = i_bucket + 1; j_bucket < bucket_size; ++j_bucket) {
        if (merged_supernode[j_bucket]) {
          continue;
        }
        const Int j_index = bucket[j_bucket];
        const Int j = supernodal_pivot_structure[j_index];
        if (graph->StructuralSupervariablesAreQuotientIndistinguishable(i, j)) {
          #pragma omp critical
          variable_merges.emplace_back(i, j, graph->supernodes[j].size());

          temp = graph->supernodes[i];
          MergeSets(temp, graph->supernodes[j], &graph->supernodes[i]);
          SwapClearVector(&graph->supernodes[j]);
          SwapClearVector(&graph->adjacency_lists[j]);
          SwapClearVector(&graph->element_lists[j]);
          merged_supernode[j_bucket] = true;
        }
      }
    }
  }

  // Update the external degree heap in a batched manner.
  graph->external_degree_heap.ReserveValueChanges(2 * variable_merges.size());
  for (const VariableMergeInfo& merge : variable_merges) {
    graph->external_degree_heap.QueueValueUpdate(
        merge.primary_index, -merge.absorbed_size);
    graph->external_degree_heap.QueueIndexDisablement(merge.absorbed_index);
  }
  graph->external_degree_heap.FlushValueChangeQueue();

  // If requested, keep track of the variable merges.
  if (store_variable_merges) {
    for (const VariableMergeInfo& merge : variable_merges) {
      graph->variable_merges.emplace_back(
          merge.primary_index, merge.absorbed_index);
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
          control.hash_type, &quotient_graph);
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
