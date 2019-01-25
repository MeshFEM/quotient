/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_MINIMUM_DEGREE_IMPL_H_
#define QUOTIENT_MINIMUM_DEGREE_IMPL_H_

#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "quotient/coordinate_graph.hpp"
#include "quotient/index_utils.hpp"
#include "quotient/integers.hpp"
#include "quotient/macros.hpp"
#include "quotient/quotient_graph.hpp"
#include "quotient/timer.hpp"

#include "quotient/minimum_degree.hpp"

namespace quotient {

inline MinimumDegreeResult::MinimumDegreeResult() {}

inline Int MinimumDegreeResult::NumStrictlyLowerCholeskyNonzeros() const {
  return num_cholesky_nonzeros - permuted_supernode_sizes.Size();
}

inline Int MinimumDegreeResult::LargestSupernode() const {
  Int largest_supernode = -1;
  Int largest_supernode_size = 0;
  for (Int i = 0; i < permuted_supernode_sizes.Size(); ++i) {
    if (permuted_supernode_sizes[i] > largest_supernode_size) {
      largest_supernode = i;
      largest_supernode_size = permuted_supernode_sizes[i];
    }
  }
  return largest_supernode;
}

inline Int MinimumDegreeResult::LargestSupernodeSize() const {
  return permuted_supernode_sizes[LargestSupernode()];
}

inline double MinimumDegreeResult::FractionOfPivotsWithMultipleElements()
    const {
  Int num_pivots_with_multiple_elements = 0;
  for (const Int& pivot_element_list_size : pivot_element_list_sizes) {
    if (pivot_element_list_size > 2) {
      ++num_pivots_with_multiple_elements;
    }
  }
  return num_pivots_with_multiple_elements /
         (1. * pivot_element_list_sizes.size());
}

inline double MinimumDegreeResult::FractionOfDegreeUpdatesWithMultipleElements()
    const {
  return num_degree_updates_with_multiple_elements / (1. * num_degree_updates);
}

inline Buffer<Int> MinimumDegreeResult::Permutation() const {
  return permutation;
}

inline Buffer<Int> MinimumDegreeResult::InversePermutation() const {
  return inverse_permutation;
}

inline void MinimumDegreeResult::PermutedAssemblyForestToDot(
    const std::string& filename) const {
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Could not open " << filename << std::endl;
    return;
  }

  file << "digraph g{\n";
  for (Int i = 0; i < permuted_assembly_parents.Size(); ++i) {
    if (permuted_assembly_parents[i] < 0) {
      continue;
    }
    std::ostringstream os;
    os << "  " << permuted_assembly_parents[i] << " -> " << i << ";\n";
    file << os.str();
  }
  file << "}\n";
}

inline MinimumDegreeResult MinimumDegree(const CoordinateGraph& graph,
                                         const MinimumDegreeControl& control) {
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

  // Eliminate the variables.
  QuotientGraph quotient_graph(graph, control);
  while (quotient_graph.NumEliminatedVertices() < num_orig_vertices) {
    quotient_graph.FindAndProcessPivot();

    if (control.store_pivot_element_list_sizes) {
      analysis.pivot_element_list_sizes.push_back(
          quotient_graph.NumPivotElements());
    }
    analysis.num_cholesky_nonzeros += quotient_graph.NumPivotCholeskyNonzeros();
    analysis.num_cholesky_flops += quotient_graph.NumPivotCholeskyFlops();
    analysis.num_degree_updates += quotient_graph.NumPivotDegreeUpdates();
    if (control.store_num_degree_updates_with_multiple_elements) {
      analysis.num_degree_updates_with_multiple_elements +=
          quotient_graph.NumPivotDegreeUpdatesWithMultipleElements();
    }
  }
  quotient_graph.CombineDenseNodes();

  // Assume the Schur complement of the non-dense supernodes onto the "dense"
  // ones results in a dense Schur complement.
  const Int num_dense = quotient_graph.NumDense();
  analysis.num_cholesky_nonzeros += ((num_dense + 1) * num_dense) / 2;
  analysis.num_cholesky_flops += std::pow(1. * num_dense, 3.) / 3.;

  // Compute the permutation using the post-ordering.
  quotient_graph.ComputePostorder(&analysis.inverse_permutation);
  InvertPermutation(analysis.inverse_permutation, &analysis.permutation);

  // Compute a map from the permuted indices to the containing supernode.
  quotient_graph.PermutedSupernodeSizes(analysis.inverse_permutation,
                                        &analysis.permuted_supernode_sizes);
  quotient_graph.PermutedMemberToSupernode(
      analysis.inverse_permutation, &analysis.permuted_member_to_supernode);

  quotient_graph.PermutedAssemblyParents(analysis.permutation,
                                         analysis.permuted_member_to_supernode,
                                         &analysis.permuted_assembly_parents);

  analysis.num_hash_collisions = quotient_graph.NumHashCollisions();
  analysis.num_hash_bucket_collisions =
      quotient_graph.NumHashBucketCollisions();
  analysis.num_aggressive_absorptions =
      quotient_graph.NumAggressiveAbsorptions();
#ifdef QUOTIENT_ENABLE_TIMERS
  const Buffer<std::pair<std::string, double>>& timings =
      quotient_graph.ComponentSeconds();
  for (const std::pair<std::string, double>& pairing : timings) {
    analysis.elapsed_seconds[pairing.first] = pairing.second;
  }
#endif

  // Extract the elimination order.
  analysis.elimination_order = quotient_graph.EliminationOrder();

  return analysis;
}

}  // namespace quotient

#endif  // ifndef QUOTIENT_MINIMUM_DEGREE_IMPL_H_
