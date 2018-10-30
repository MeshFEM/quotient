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
#include "quotient/macros.hpp"
#include "quotient/coordinate_graph.hpp"
#include "quotient/minimum_degree.hpp"
#include "quotient/quotient_graph.hpp"
#include "quotient/timer.hpp"

namespace quotient {

inline MinimumDegreeResult::MinimumDegreeResult() { }

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

inline std::vector<Int>
MinimumDegreeResult::Permutation() const {
  const Int num_vertices = preorder.size();
#ifdef QUOTIENT_DEBUG
  std::vector<Int> permutation(num_vertices, -1);
#else
  std::vector<Int> permutation(num_vertices);
#endif

  // Fill the permutation with the inverse of the postordering, which is the
  // reverse of the preordering.
  for (Int index = 0; index < num_vertices; ++index) {
    permutation[preorder[num_vertices - 1 - index]] = index;
  }

#ifdef QUOTIENT_DEBUG
  for (Int index = 0; index < num_vertices; ++index) {
    QUOTIENT_ASSERT(permutation[index] != -1,
        "Permutation was only partially filled.");
  }
#endif

  return permutation;
}

inline MinimumDegreeResult MinimumDegree(
  const CoordinateGraph& graph, const MinimumDegreeControl& control) {
  QUOTIENT_ASSERT(graph.NumSources() == graph.NumTargets(),
      "MinimumDegree requires a symmetric input graph.");
  const Int num_orig_vertices = graph.NumSources();

  // Initialize a data structure that will eventually contain the results of
  // the (approximae) minimum degree analysis.
  MinimumDegreeResult analysis;
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
  external_degrees.reserve(num_orig_vertices - 1);

  // A vector for storing the hashes of the supervariables in the current
  // pivot's structure.
  std::vector<std::size_t> bucket_keys;
  bucket_keys.reserve(num_orig_vertices - 1);

  // Set up a set of timers for the components of the analysis.
  std::unordered_map<std::string, Timer> timers;
  constexpr char kSetup[] = "Setup";
  constexpr char kComputePivotStructure[] = "ComputePivotStructure";
  constexpr char kAbsorption[] = "Absorption";
  constexpr char kResetExternalElementSizes[] = "ResetExternalElementSizes";
  constexpr char kComputeExternalDegrees[] = "ComputeExternalDegrees";
  constexpr char kUpdateExternalDegrees[] = "UpdateExternalDegrees";
  constexpr char kMergeVariables[] = "MergeVariables";
  constexpr char kFinalize[] = "Finalize";

  // Eliminate the variables.
  if (control.time_stages) timers[kSetup].Start();
  QuotientGraph quotient_graph(graph, control);
  if (control.time_stages) timers[kSetup].Stop();
  while (quotient_graph.NumEliminatedVertices() < num_orig_vertices) {
    // Get the next pivot.
    quotient_graph.GetNextPivot();
    if (control.store_pivot_element_list_sizes) {
      analysis.pivot_element_list_sizes.push_back(
          quotient_graph.NumPivotElements());
    }

    // Compute the structure of this pivot block.
    if (control.time_stages) timers[kComputePivotStructure].Start();
    quotient_graph.ComputePivotStructure();
    if (control.time_stages) timers[kComputePivotStructure].Stop();
    analysis.num_cholesky_nonzeros += quotient_graph.NumPivotCholeskyNonzeros();
    analysis.num_cholesky_flops += quotient_graph.NumPivotCholeskyFlops();

    // Compute the external structure cardinalities, |L_e \ L_p|, of all
    // elements e in an element list of a supernode in L_p. Any elements to
    // be aggressively absorbed are also returned.
    if (control.time_stages) timers[kAbsorption].Start();
    quotient_graph.AbsorptionAndExternalElementSizes(
        &aggressive_absorption_elements);
    if (control.time_stages) timers[kAbsorption].Stop();

    // Store the external degrees of all supervariables in the pivot structure.
    if (control.time_stages) timers[kComputeExternalDegrees].Start();
    quotient_graph.ComputeExternalDegreesAndHashes(
        &external_degrees, &bucket_keys);
    if (control.time_stages) timers[kComputeExternalDegrees].Stop();

    // Update the degree lists using the computed degrees.
    if (control.time_stages) timers[kUpdateExternalDegrees].Start();
    quotient_graph.UpdateExternalDegrees(external_degrees);
    if (control.time_stages) timers[kUpdateExternalDegrees].Stop();

    // Update/store metadata associated with the degree computations.
    analysis.num_degree_updates += quotient_graph.NumPivotDegreeUpdates();
    if (control.store_num_degree_updates_with_multiple_elements) {
      analysis.num_degree_updates_with_multiple_elements +=
          quotient_graph.NumPivotDegreeUpdatesWithMultipleElements();
    }

    if (control.allow_supernodes) {
      // Merge any equivalent supernodes by explicitly checking for equality
      // between pairs that are in the same hash bucket.
      if (control.time_stages) timers[kMergeVariables].Start();
      quotient_graph.MergeVariables(bucket_keys);
      if (control.time_stages) timers[kMergeVariables].Stop();
    }

    // Clear the external element size array.
    if (control.time_stages) timers[kResetExternalElementSizes].Start();
    quotient_graph.ResetExternalElementSizes(); 
    if (control.time_stages) timers[kResetExternalElementSizes].Stop();

    // Formally convert the pivot from a supervariable into an element.
    quotient_graph.ConvertPivotIntoElement(aggressive_absorption_elements);
  }

  // Extract the relevant information from the QuotientGraph.
  if (control.time_stages) timers[kFinalize].Start();
  analysis.supernodes.resize(num_orig_vertices);
  for (Int i = 0; i < num_orig_vertices; ++i) {
    analysis.supernodes[i] = quotient_graph.FormSupernode(i);
  }
  analysis.elimination_order = quotient_graph.EliminationOrder();
  quotient_graph.ComputePreorder(&analysis.preorder);
  if (control.store_structures) {
    quotient_graph.FormEliminatedStructures(&analysis.eliminated_structures);
  }
  analysis.num_hash_collisions = quotient_graph.NumHashCollisions();
  analysis.num_hash_bucket_collisions =
      quotient_graph.NumHashBucketCollisions();
  analysis.num_aggressive_absorptions =
      quotient_graph.NumAggressiveAbsorptions();
  for (const std::pair<std::string, Timer>& pairing : timers) {
    analysis.elapsed_seconds[pairing.first] = pairing.second.TotalSeconds();
  }
  if (control.time_stages) timers[kFinalize].Stop();

  return analysis;
}

} // namespace quotient

#endif // ifndef QUOTIENT_MINIMUM_DEGREE_IMPL_H_
