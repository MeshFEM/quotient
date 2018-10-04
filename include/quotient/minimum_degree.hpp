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

namespace quotient {

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

  // The supernodal nonzero structure of the principal columns of the
  // factorization.
  std::vector<std::vector<Int>> supernodal_structures;

  // We will push to elimination order as the reordering algorithm progresses,
  // so we will allocate an upper bound for the amount of required space.
  // The 'supernodes' and 'structures' variables will be copied over from
  // the quotient graph just before the analysis completes.
  MinimumDegreeAnalysis(Int num_vertices) {
    elimination_order.reserve(num_vertices);
  }
};

// Compute the structure of the pivot:
//   L_p := (A_p \cup (\cup_{e in E_p} L_e)) \ p.
//
// TODO(Jack Poulson): Experiment with several different union-finding
// algorithms. It isn't clear that repeated calls to std::set_union has
// optimal time complexity.
inline void ComputePivotStructure(Int pivot, QuotientGraph* graph) {
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

// Update the adjacency lists after computing the supernodal pivot structure.
inline void UpdateAdjacencyListsAfterSelectingPivot(
    Int pivot,
    const std::vector<Int>& supernodal_pivot_structure,
    QuotientGraph* graph) {
  std::vector<Int> temp;
  for (const Int& i : supernodal_pivot_structure) {
    // Remove redundant adjacency entries:
    //   A_i := (A_i \ L_p) \ supernode(p).
    FilterSet(graph->adjacency_lists[i], graph->structures[pivot], &temp);
    FilterSet(temp, graph->supernodes[pivot], &graph->adjacency_lists[i]);
  }
}

// Update the element lists after computing the supernodal pivot structure.
inline void UpdateElementListsAfterSelectingPivot(
    Int pivot,
    const std::vector<Int>& supernodal_pivot_structure,
    QuotientGraph* graph) {
  std::vector<Int> temp;
  for (const Int& i : supernodal_pivot_structure) {
    // Element absorption:
    //   E_i := (E_i \ E_p) \cup {p}.
    temp = graph->element_lists[i];
    FilterSet(temp, graph->element_lists[pivot], &graph->element_lists[i]);
    InsertEntryIntoSet(pivot, &graph->element_lists[i]);
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
    ExternalDegreeType degree_type,
    QuotientGraph* graph) {
  // Compute the external structure cardinalities of the elements.
  // (but only if the Amestoy external degree approximation is requested).
  std::unordered_map<Int, Int> external_structure_sizes;
  if (degree_type == kAmestoyExternalDegree) {
    external_structure_sizes = graph->ExternalStructureSizes(
        supernodal_pivot_structure);
  }

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
    const std::size_t bucket = graph->AshcraftVariableHash(i);
    bucket_list[i_index] = bucket;
    variable_hash_map[bucket].push_back(i_index);
  }

  std::vector<bool> merged_supernode(supernodal_pivot_struct_size, false);
  std::vector<Int> temp;
  for (Int i_index = 0; i_index < supernodal_pivot_struct_size; ++i_index) {
    if (merged_supernode[i_index]) {
      continue;
    }
    const Int i = supernodal_pivot_structure[i_index];
    const Int bucket = bucket_list[i_index];
    for (const Int& j_index : variable_hash_map[bucket]) {
      if (j_index <= i_index || merged_supernode[j_index]) {
        continue;
      }
      const Int j = supernodal_pivot_structure[j_index];
      if (graph->VariablesAreQuotientIndistinguishable(i, j)) {
        // Absorb supernode(j) into supernode(i). 
        temp = graph->supernodes[i];
        MergeSets(temp, graph->supernodes[j], &graph->supernodes[i]);
        graph->external_degree_heap.UpdateValue(
            i, -graph->supernodes[j].size());
        graph->external_degree_heap.DisableIndex(j);
#ifdef QUOTIENT_DEBUG
        graph->variables.erase(j);
#endif
        SwapClearVector(&graph->supernodes[j]);
        SwapClearVector(&graph->adjacency_lists[j]);
        SwapClearVector(&graph->element_lists[j]);
        merged_supernode[j_index] = true;
      }
    }
  }
}

// Converts the 'pivot' (super)variable into an element.
inline void ConvertPivotIntoElement(Int pivot, QuotientGraph* graph) {
#ifdef QUOTIENT_DEBUG
  // Update the element list:
  //   \bar{V} := (\bar{V} \cup {p}) \ E_p
  graph->elements.insert(pivot);
  for (const Int& element : graph->element_lists[pivot]) {
    graph->elements.erase(element);
  }
  // Update the variable list:
  //   V := V \ {p}.
  graph->variables.erase(pivot);
#endif
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
  const CoordinateGraph& graph, ExternalDegreeType degree_type) {
#ifdef QUOTIENT_DEBUG
  if (graph.NumSources() != graph.NumVertices()) {
    std::cerr << "ERROR: MinimumDegree requires a symmetric input graph."
              << std::endl;
    return MinimumDegreeAnalysis(0);
  }
#endif
  QuotientGraph quotient_graph(graph);
  const Int num_orig_vertices = quotient_graph.num_original_vertices;

  // Eliminate the variables.
  MinimumDegreeAnalysis analysis(num_orig_vertices);
  while (quotient_graph.num_eliminated_vertices < num_orig_vertices) {
    // Retrieve a variable with minimal (approximate) external degree.
    const std::pair<Int, Int> pivot_pair =
        quotient_graph.external_degree_heap.MinimalEntry(); 
    const Int pivot = pivot_pair.first;

    ComputePivotStructure(pivot, &quotient_graph);
    const std::vector<Int> supernodal_pivot_structure =
        quotient_graph.FormSupernodalStructure(pivot);

    UpdateAdjacencyListsAfterSelectingPivot(
        pivot, supernodal_pivot_structure, &quotient_graph);

    UpdateElementListsAfterSelectingPivot(
        pivot, supernodal_pivot_structure, &quotient_graph);

    UpdateExternalDegreeApproximations(
        pivot, supernodal_pivot_structure, degree_type, &quotient_graph);

    DetectAndMergeVariables(supernodal_pivot_structure, &quotient_graph);

    ConvertPivotIntoElement(pivot, &quotient_graph);

    analysis.elimination_order.push_back(pivot);
  }

  analysis.supernodes = quotient_graph.supernodes;
  analysis.supernodal_structures.resize(num_orig_vertices);
  for (const Int& i : analysis.elimination_order) {
    analysis.supernodal_structures[i] =
        quotient_graph.FormSupernodalStructure(i);
  }
  return analysis;
}

} // namespace quotient

#endif // ifndef QUOTIENT_MINIMUM_DEGREE_H_
