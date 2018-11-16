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

#include "quotient/integers.hpp"
#include "quotient/macros.hpp"
#include "quotient/coordinate_graph.hpp"
#include "quotient/quotient_graph.hpp"
#include "quotient/timer.hpp"

#include "quotient/minimum_degree.hpp"

namespace quotient {

inline MinimumDegreeResult::MinimumDegreeResult() { }

inline Int MinimumDegreeResult::NumStrictlyLowerCholeskyNonzeros() const {
  return num_cholesky_nonzeros - supernode_sizes.size();
}

inline Int MinimumDegreeResult::LargestSupernode() const {
  Int largest_supernode = -1;
  Int largest_supernode_size = 0;
  for (std::size_t i = 0; i < supernode_sizes.size(); ++i) {
    if (supernode_sizes[i] > largest_supernode_size) {
      largest_supernode = i;
      largest_supernode_size = supernode_sizes[i];
    }
  }
  return largest_supernode;
}

inline Int MinimumDegreeResult::LargestSupernodeSize() const {
  return supernode_sizes[LargestSupernode()];
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
  const Int num_vertices = postorder.size();
#ifdef QUOTIENT_DEBUG
  std::vector<Int> permutation(num_vertices, -1);
#else
  std::vector<Int> permutation(num_vertices);
#endif

  // Fill the permutation with the inverse of the postordering.
  for (Int index = 0; index < num_vertices; ++index) {
    permutation[postorder[index]] = index;
  }

#ifdef QUOTIENT_DEBUG
  for (Int index = 0; index < num_vertices; ++index) {
    QUOTIENT_ASSERT(permutation[index] != -1,
        "Permutation was only partially filled.");
  }
#endif

  return permutation;
}

inline void MinimumDegreeResult::AssemblyForestToDot(
    const std::string& filename) const {
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Could not open " << filename << std::endl;
    return;
  }

  file << "digraph g{\n";
  for (const Int& i : postorder) {
    // Skip empty and root nodes.
    if (!supernode_sizes[i] || assembly_parents[i] <= 0) {
      continue;
    }
    std::ostringstream os;
    os << "  " << assembly_parents[i] << " -> " << i << ";\n";
    file << os.str();
  }
  file << "}\n";
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

  // Extract the relevant information from the QuotientGraph.
  if (control.store_supernodes) {
    analysis.supernodes.resize(num_orig_vertices);
    for (Int i = 0; i < num_orig_vertices; ++i) {
      analysis.supernodes[i] = quotient_graph.FormSupernode(i);
    }
  }
  analysis.supernode_sizes.resize(num_orig_vertices);
  for (Int i = 0; i < num_orig_vertices; ++i) {
    analysis.supernode_sizes[i] = quotient_graph.SupernodeSize(i);
  }
  analysis.elimination_order = quotient_graph.EliminationOrder();
  quotient_graph.ComputePostorder(&analysis.postorder);
  analysis.assembly_parents = quotient_graph.AssemblyParents();
  if (control.store_structures) {
    quotient_graph.FormEliminatedStructures(&analysis.eliminated_structures);
  }
  analysis.num_hash_collisions = quotient_graph.NumHashCollisions();
  analysis.num_hash_bucket_collisions =
      quotient_graph.NumHashBucketCollisions();
  analysis.num_aggressive_absorptions =
      quotient_graph.NumAggressiveAbsorptions();
#ifdef QUOTIENT_ENABLE_TIMERS
  const std::vector<std::pair<std::string, double>>& timings =
      quotient_graph.ComponentSeconds();
  for (const std::pair<std::string, double>& pairing : timings) {
    analysis.elapsed_seconds[pairing.first] = pairing.second;
  }
#endif

  return analysis;
}

} // namespace quotient

#endif // ifndef QUOTIENT_MINIMUM_DEGREE_IMPL_H_
