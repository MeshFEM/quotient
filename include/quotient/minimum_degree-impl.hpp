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

// Computes the structure of the pivot:
//
//   L_p := (A_p \cup (\cup_{e in E_p} L_e)) \ supernode(p).
//
// It is assumed that the mask is of length 'num_orig_vertices' and set to all
// zeros on input.
inline void ComputePivotStructure(
    Int pivot,
    QuotientGraph* graph,
    std::vector<Int>* supernodal_pivot_structure,
    std::vector<int>* pivot_mask) {
  graph->structures[pivot].clear();
  supernodal_pivot_structure->clear();
  for (const Int& index : graph->adjacency_lists[pivot]) {
    if (graph->supernode_sizes[index] == 0) {
      continue;
    }
    (*pivot_mask)[index] = 1;
    supernodal_pivot_structure->push_back(index);

    // Fill in the supernode.
    Int k = index;
    const Int k_last = graph->tail_index[index];
    while (true) {
      graph->structures[pivot].push_back(k);
      if (k == k_last) {
        break;
      }
      k = graph->next_index[k];
    }
  }
  for (const Int& element : graph->element_lists[pivot]) {
    for (const Int& index : graph->structures[element]) {
      if (index == pivot ||
          graph->supernode_sizes[index] == 0 ||
          (*pivot_mask)[index]) {
        continue;
      }
      (*pivot_mask)[index] = 1;
      supernodal_pivot_structure->push_back(index);

      // Fill in the supernode.
      Int k = index;
      const Int k_last = graph->tail_index[index];
      while (true) {
        graph->structures[pivot].push_back(k);
        if (k == k_last) {
          break;
        }
        k = graph->next_index[k];
      }
    }
  }
}

// Sets the entry mask[i] to zero for each i in indices.
inline void ClearMask(const std::vector<Int>& indices, std::vector<int>* mask) {
  for (const Int& i : indices) {
    (*mask)[i] = 0;
  }
}

inline bool MaskIsClear(const std::vector<int>& mask) {
  for (const int& value : mask) {
    if (value) {
      return false;
    }
  }
  return true;
}

// Update the adjacency lists after computing the supernodal pivot structure.
inline void UpdateAdjacencyListsAfterSelectingPivot(
    Int pivot,
    const std::vector<Int>& supernodal_pivot_structure,
    const std::vector<int>& pivot_mask,
    QuotientGraph* graph) {
  for (const Int& i : supernodal_pivot_structure) {
    // Remove redundant adjacency entries:
    //   A_i := (A_i \ L_p) \ supernode(p).
    Int packed_size = 0;
    for (Int index : graph->adjacency_lists[i]) {
      if (index == pivot ||
          pivot_mask[index] ||
          graph->supernode_sizes[index] == 0) {
        continue;
      }
      graph->adjacency_lists[i][packed_size++] = index;
    }
    graph->adjacency_lists[i].resize(packed_size);
  }
}

// Update the element lists after computing the supernodal pivot structure.
inline void UpdateElementListsAfterSelectingPivot(
    Int pivot,
    const std::vector<Int>& supernodal_pivot_structure,
    const std::vector<Int>& aggressive_absorption_elements,
    bool store_aggressive_absorptions,
    QuotientGraph* graph,
    std::vector<int>* pivot_mask) {
  // Mark the pivot elements in the mask.
  for (const Int& element : graph->element_lists[pivot]) {
    (*pivot_mask)[element] = -1;
  }

  for (const Int& i : supernodal_pivot_structure) {
    // Element absorption:
    //   E_i := (E_i \ E_p) \cup {p}.
    std::vector<Int>& element_list = graph->element_lists[i];
    Int num_packed = 0;
    for (const Int& element : element_list) {
      if ((*pivot_mask)[element]) { 
        continue;
      }
      element_list[num_packed++] = element;
    }
    element_list.resize(num_packed);
    element_list.push_back(pivot);
  }

  if (!aggressive_absorption_elements.empty()) {
    for (const Int& element : aggressive_absorption_elements) {
      // Aggressive element absorption:
      //   E_p := (E_p \ E_e) \cup {e}.
      // Instead of explicitly forming this update, we will disable the
      // elements from E_e in the pivot_mask and implicitly keep track of
      // the fact that {e} should be inserted.
      if (store_aggressive_absorptions) {
        graph->aggressive_absorptions.emplace_back(pivot, element);
      }
      for (const Int& e : graph->element_lists[element]) {
        (*pivot_mask)[e] = 0;
      }
    }

    // Form E_p.
    Int num_packed = 0;
    std::vector<Int>& element_list = graph->element_lists[pivot];
    for (const Int& element : element_list) {
      if ((*pivot_mask)[element]) {
        element_list[num_packed++] = element;
      }
    }
    element_list.resize(num_packed);
    for (const Int& element : aggressive_absorption_elements) {
      if ((*pivot_mask)[element]) {
        continue;
      }
      element_list.push_back(element);
    }
  }

  // Unset the pivot elements in the mask.
  for (const Int& element : graph->element_lists[pivot]) {
    (*pivot_mask)[element] = 0;
  }
}

// Compute the external degree approximations of the supernodes
// adjacent to the current pivot.
inline void UpdateExternalDegreeApproximations(
    Int pivot,
    const std::vector<Int>& supernodal_pivot_structure,
    const std::vector<int>& pivot_mask,
    const std::vector<Int>& external_structure_sizes,
    ExternalDegreeType degree_type,
    std::vector<int>* exact_degree_mask,
    QuotientGraph* graph) {
  // Compute and store the external degrees.
  const std::size_t supernodal_struct_size = supernodal_pivot_structure.size();
  std::vector<Int> external_degrees(supernodal_struct_size);
  // TODO(Jack Poulson): Consider how to parallelize using different masks
  // for each thread.
  for (std::size_t index = 0; index < supernodal_struct_size; ++index) {
    const Int i = supernodal_pivot_structure[index];

    // Compute the external degree (or approximation) of supervariable i:
    //   d_i := |A_i \ supernode(i)| +
    //       |(\cup_{e in E_i} L_e) \ supernode(i)|.
    external_degrees[index] = ExternalDegree(
        *graph, i, pivot, pivot_mask, external_structure_sizes, degree_type,
        exact_degree_mask);
  }

  // Insert the external degrees into the heap.
  graph->external_degree_heap.ReserveValueChanges(supernodal_struct_size);
  for (std::size_t index = 0; index < supernodal_struct_size; ++index) {
    graph->external_degree_heap.QueueValueAssignment(
        supernodal_pivot_structure[index], external_degrees[index]);
  }
  graph->external_degree_heap.FlushValueChangeQueue();
}

// Computes hashes of the supervariables in the pivot structure.
inline std::vector<std::size_t> ComputeVariableHashes(
    const QuotientGraph& graph,
    const std::vector<Int>& supernodal_pivot_structure,
    VariableHashType hash_type) {
  // Compute the hashes for each variable.
  const std::size_t supernodal_struct_size = supernodal_pivot_structure.size();
  std::vector<std::size_t> bucket_keys(supernodal_struct_size);
  #pragma omp parallel for schedule(dynamic)
  for (std::size_t i_index = 0; i_index < supernodal_struct_size; ++i_index) {
    const Int i = supernodal_pivot_structure[i_index];
    bucket_keys[i_index] = graph.VariableHash(i, hash_type);
  }
  return bucket_keys;
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
    const std::vector<std::size_t>& bucket_keys,
    bool store_variable_merges,
    QuotientGraph* graph,
    std::vector<std::vector<Int>>* buckets) {
  // Fill a set of buckets for the hashes of the supernodes adjacent to
  // the current pivot.
  const Int supernodal_struct_size = bucket_keys.size();
  std::vector<Int> bucket_indices;
  for (Int i_index = 0; i_index < supernodal_struct_size; ++i_index) {
    const std::size_t bucket_key = bucket_keys[i_index];
    const Int bucket_index = bucket_key % graph->num_original_vertices;
    (*buckets)[bucket_index].push_back(i_index);
    if ((*buckets)[bucket_index].size() == 2) {
      bucket_indices.push_back(bucket_index);
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

  std::vector<VariableMergeInfo> variable_merges;
  #pragma omp parallel for schedule(dynamic)
  for (std::size_t index = 0; index < bucket_indices.size(); ++index) {
    const Int bucket_index = bucket_indices[index];
    const std::vector<Int>& bucket = (*buckets)[bucket_index];
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
          variable_merges.emplace_back(i, j, graph->supernode_sizes[j]);

          // Merge [i] -> [j] (and i becomes the principal member).
          Int k = j;
          while (true) {
            graph->head_index[k] = i;
            if (k == graph->tail_index[j]) {
              break;
            }
            k = graph->next_index[k];
          }
          graph->next_index[graph->tail_index[i]] = j;
          graph->tail_index[i] = graph->tail_index[j];
          graph->supernode_sizes[i] += graph->supernode_sizes[j];
          graph->supernode_sizes[j] = 0;

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

  // Clear the buckets.
  for (Int i_index = 0; i_index < supernodal_struct_size; ++i_index) {
    const std::size_t bucket_key = bucket_keys[i_index];
    const Int bucket_index = bucket_key % graph->num_original_vertices;
    (*buckets)[bucket_index].clear();
  }
}

// Converts the 'pivot' (super)variable into an element.
inline void ConvertPivotIntoElement(Int pivot, QuotientGraph* graph) {
  graph->external_degree_heap.DisableIndex(pivot);
  SwapClearVector(&graph->adjacency_lists[pivot]);
  SwapClearVector(&graph->element_lists[pivot]);
  graph->num_eliminated_vertices += graph->supernode_sizes[pivot];
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
  const Int num_orig_vertices = graph.NumSources();

  // Initialize a data structure that will eventually contain the results of
  // the (approximae) minimum degree analysis.
  MinimumDegreeAnalysis analysis(num_orig_vertices);
  if (control.store_pivot_element_list_sizes) {
    analysis.pivot_element_list_sizes.reserve(num_orig_vertices);
  }
  if (control.store_num_degree_updates_with_multiple_elements) {
    analysis.num_degree_updates_with_multiple_elements = 0;
    analysis.num_degree_updates_without_multiple_elements = 0;
  }

  // A mask of length 'num_orig_vertices' that can be used to quickly compute
  // the cardinalities of |L_e \ L_p| for each element e in an element list of
  // a supervariable in the current pivot structure, L_p.
  std::vector<Int> external_structure_sizes;
  const bool compute_external_structure_sizes =
      control.aggressive_absorption ||
      control.degree_type == kAmestoyExternalDegree;

  // A vector that will be used to store the list of elements that should be
  // aggressively absorbed in a particular stage.
  std::vector<Int> aggressive_absorption_elements;

  // A list of the principal members of the current pivot's structure.
  std::vector<Int> supernodal_pivot_structure;

  // A mask of length 'num_orig_vertices' that is 1 in index 'i' if and only
  // if 'i' is a member of the current pivot's structure. All other entries
  // will be zero.
  std::vector<int> pivot_mask(num_orig_vertices, 0);

  // A mask of length 'num_orig_vertices' that used within exact external
  // degree computations to perform set unions. It is only created if exact
  // degree computations were requested, and it must be set to all zeros before
  // and after each call to ExternalDegree.
  std::vector<int> exact_degree_mask;
  if (control.degree_type == kExactExternalDegree) {
    exact_degree_mask.resize(num_orig_vertices, 0);
  }

  // Set up a set of timers for the components of the analysis.
  std::unordered_map<std::string, Timer> timers;
  constexpr char kComputePivotStructure[] = "ComputePivotStructure";
  constexpr char kUpdateAdjacencyLists[] = "UpdateAdjacencyLists";
  constexpr char kExternalStructureSizes[] = "ExternalStructureSizes";
  constexpr char kUpdateElementLists[] = "UpdateElementLists";
  constexpr char kUpdateExternalDegrees[] = "UpdateExternalDegrees";
  constexpr char kComputeVariableHashes[] = "ComputeVariableHashes";
  constexpr char kDetectAndMergeVariables[] = "DetectAndMergeVariables";
  if (control.time_stages) {
    timers[kComputePivotStructure].Reset();
    timers[kUpdateAdjacencyLists].Reset();
    timers[kExternalStructureSizes].Reset();
    timers[kUpdateElementLists].Reset();
    timers[kUpdateExternalDegrees].Reset();
    timers[kComputeVariableHashes].Reset();
    timers[kDetectAndMergeVariables].Reset();
  }

  // A set of buckets for each hash value (modulo num_original_vertices) of the
  // supervariables.
  std::vector<std::vector<Int>> buckets(num_orig_vertices);

  // Eliminate the variables.
  QuotientGraph quotient_graph(graph);
  if (compute_external_structure_sizes) {
    quotient_graph.InitializeExternalStructureSizes(&external_structure_sizes);
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
    ComputePivotStructure(
        pivot, &quotient_graph, &supernodal_pivot_structure, &pivot_mask);
    if (control.time_stages) timers[kComputePivotStructure].Stop();

    if (control.time_stages) timers[kUpdateAdjacencyLists].Start();
    UpdateAdjacencyListsAfterSelectingPivot(
        pivot, supernodal_pivot_structure, pivot_mask, &quotient_graph);
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
        control.store_aggressive_absorptions, &quotient_graph, &pivot_mask);
    if (control.time_stages) timers[kUpdateElementLists].Stop();

    if (control.time_stages) timers[kUpdateExternalDegrees].Start();
    UpdateExternalDegreeApproximations(
        pivot, supernodal_pivot_structure, pivot_mask,
        external_structure_sizes, control.degree_type, &exact_degree_mask,
        &quotient_graph);
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
    ClearMask(supernodal_pivot_structure, &pivot_mask);

    if (control.allow_supernodes) {
      if (control.time_stages) timers[kComputeVariableHashes].Start();
      const std::vector<std::size_t> bucket_keys = ComputeVariableHashes(
          quotient_graph, supernodal_pivot_structure, control.hash_type);
      if (control.time_stages) timers[kComputeVariableHashes].Stop();
      if (control.time_stages) timers[kDetectAndMergeVariables].Start();
      DetectAndMergeVariables(
          supernodal_pivot_structure, bucket_keys,
          control.store_variable_merges, &quotient_graph, &buckets);
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
  analysis.supernodes.resize(num_orig_vertices);
  for (Int i = 0; i < num_orig_vertices; ++i) {
    if (quotient_graph.supernode_sizes[i] > 0) {
      analysis.supernodes[i] = quotient_graph.FormSupernode(i);
    }
  }
  analysis.principal_structures.resize(analysis.elimination_order.size());
  for (std::size_t index = 0; index < analysis.elimination_order.size();
       ++index) {
    analysis.principal_structures[index] = 
        quotient_graph.structures[analysis.elimination_order[index]];
    std::sort(
        analysis.principal_structures[index].begin(),
        analysis.principal_structures[index].end());
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
