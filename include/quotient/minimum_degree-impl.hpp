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
#include "quotient/minimum_degree.hpp"
#include "quotient/quotient_graph.hpp"
#include "quotient/timer.hpp"

namespace quotient {

inline MinimumDegreeResult::MinimumDegreeResult(Int num_vertices) {
  elimination_order.reserve(num_vertices);
}

inline Int MinimumDegreeResult::NumStrictlyLowerCholeskyNonzeros() const {
  return num_cholesky_nonzeros - supernodes.size();
}

inline Int MinimumDegreeResult::LargestSupernode() const {
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

inline Int MinimumDegreeResult::LargestSupernodeSize() const {
  return supernodes[LargestSupernode()].size();
}

inline double
MinimumDegreeResult::FractionOfPivotsWithMultipleElements() const {
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
MinimumDegreeResult::FractionOfDegreeUpdatesWithMultipleElements() const {
  return num_degree_updates_with_multiple_elements / (1. * num_degree_updates);
}

inline MinimumDegreeResult MinimumDegree(
  const CoordinateGraph& graph, const MinimumDegreeControl& control) {
#ifdef QUOTIENT_DEBUG
  if (graph.NumSources() != graph.NumTargets()) {
    std::cerr << "ERROR: MinimumDegree requires a symmetric input graph."
              << std::endl;
    return MinimumDegreeResult(0);
  }
#endif
  const Int num_orig_vertices = graph.NumSources();

  // Initialize a data structure that will eventually contain the results of
  // the (approximae) minimum degree analysis.
  MinimumDegreeResult analysis(num_orig_vertices);
  if (control.store_pivot_element_list_sizes) {
    analysis.pivot_element_list_sizes.reserve(num_orig_vertices);
  }
  if (control.store_num_degree_updates_with_multiple_elements) {
    analysis.num_degree_updates_with_multiple_elements = 0;
  }

  // A vector that will be used to store the list of elements that should be
  // aggressively absorbed in a particular stage.
  std::vector<Int> aggressive_absorption_elements;
  if (control.aggressive_absorption) {
    aggressive_absorption_elements.reserve(num_orig_vertices - 1);
  }

  // A vector for storing the list of new external degree updates.
  std::vector<Int> external_degrees;

  // A vector for storing the hashes of the supervariables in the current
  // pivot's structure.
  std::vector<std::size_t> bucket_keys;

  // Set up a set of timers for the components of the analysis.
  std::unordered_map<std::string, Timer> timers;
  constexpr char kComputePivotStructure[] = "ComputePivotStructure";
  constexpr char kUpdateAdjacencyLists[] = "UpdateAdjacencyLists";
  constexpr char kExternalElementSizes[] = "ExternalElementSizes";
  constexpr char kResetExternalElementSizes[] = "ResetExternalElementSizes";
  constexpr char kUpdateElementLists[] = "UpdateElementLists";
  constexpr char kComputeExternalDegrees[] = "ComputeExternalDegrees";
  constexpr char kUpdateExternalDegrees[] = "UpdateExternalDegrees";
  constexpr char kComputeVariableHashes[] = "ComputeVariableHashes";
  constexpr char kMergeVariables[] = "MergeVariables";

  // Eliminate the variables.
  QuotientGraph quotient_graph(graph, control);
  while (quotient_graph.NumEliminatedVertices() < num_orig_vertices) {
    const Int pivot = quotient_graph.GetNextPivot();
    if (control.store_pivot_element_list_sizes) {
      analysis.pivot_element_list_sizes.push_back(
          quotient_graph.NumPivotElements());
    }

    if (control.time_stages) timers[kComputePivotStructure].Start();
    const Int num_stale_element_members =
        quotient_graph.ComputePivotStructure();
    if (control.time_stages) timers[kComputePivotStructure].Stop();
    analysis.num_stale_element_members += num_stale_element_members;
    analysis.num_cholesky_nonzeros += quotient_graph.NumPivotCholeskyNonzeros();
    analysis.num_cholesky_flops += quotient_graph.NumPivotCholeskyFlops();

    if (control.time_stages) timers[kUpdateAdjacencyLists].Start();
    quotient_graph.UpdateAdjacencyListsAfterSelectingPivot();
    if (control.time_stages) timers[kUpdateAdjacencyLists].Stop();

    if (control.time_stages) timers[kUpdateElementLists].Start();
    quotient_graph.FlagPivotElementList();
    quotient_graph.UpdateElementListsAfterSelectingPivot();
    if (control.time_stages) timers[kUpdateElementLists].Stop();

    // Compute the external structure cardinalities of the elements.
    // (but only if the Amestoy external degree approximation is requested).
    if (quotient_graph.UsingExternalElementSizes()) {
      if (control.time_stages) timers[kExternalElementSizes].Start();
      quotient_graph.RecomputeExternalElementSizes(
          &aggressive_absorption_elements);
      if (control.time_stages) timers[kExternalElementSizes].Stop();
    }

    if (control.time_stages) timers[kUpdateElementLists].Start();
    quotient_graph.AggressiveAbsorption(aggressive_absorption_elements);
    if (control.time_stages) timers[kUpdateElementLists].Stop();

    quotient_graph.UnflagPivotElementList();

    if (control.time_stages) timers[kComputeExternalDegrees].Start();
    quotient_graph.ComputeExternalDegrees(&external_degrees);
    if (control.time_stages) timers[kComputeExternalDegrees].Stop();

    if (control.time_stages) timers[kUpdateExternalDegrees].Start();
    quotient_graph.UpdateExternalDegrees(external_degrees);
    if (control.time_stages) timers[kUpdateExternalDegrees].Stop();

    analysis.num_degree_updates += quotient_graph.NumPivotDegreeUpdates();
    if (control.store_num_degree_updates_with_multiple_elements) {
      analysis.num_degree_updates_with_multiple_elements +=
          quotient_graph.NumPivotDegreeUpdatesWithMultipleElements();
    }

    quotient_graph.UnflagPivotStructure();

    if (control.allow_supernodes) {
      if (control.time_stages) timers[kComputeVariableHashes].Start();
      quotient_graph.ComputeVariableHashes(&bucket_keys);
      if (control.time_stages) timers[kComputeVariableHashes].Stop();

      if (control.time_stages) timers[kMergeVariables].Start();
      quotient_graph.MergeVariables(bucket_keys);
      if (control.time_stages) timers[kMergeVariables].Stop();
    }

    if (quotient_graph.UsingExternalElementSizes()) {
      if (control.time_stages) timers[kResetExternalElementSizes].Start();
      quotient_graph.ResetExternalElementSizes(); 
      if (control.time_stages) timers[kResetExternalElementSizes].Stop();
    }

    quotient_graph.ConvertPivotIntoElement();

    analysis.elimination_order.push_back(pivot);
  }

  // Extract the relevant information from the QuotientGraph.
  analysis.supernodes.resize(num_orig_vertices);
  for (Int i = 0; i < num_orig_vertices; ++i) {
    analysis.supernodes[i] = quotient_graph.FormSupernode(i);
  }
  if (control.store_structures) {
    quotient_graph.FormEliminatedStructures(&analysis.eliminated_structures);
  }
  analysis.num_hash_collisions = quotient_graph.NumHashCollisions();
  analysis.aggressive_absorptions = quotient_graph.AggressiveAbsorptions();
  analysis.variable_merges = quotient_graph.VariableMerges();
  for (const std::pair<std::string, Timer>& pairing : timers) {
    analysis.elapsed_seconds[pairing.first] = pairing.second.TotalSeconds();
  }

  return analysis;
}

} // namespace quotient

#endif // ifndef QUOTIENT_MINIMUM_DEGREE_IMPL_H_
