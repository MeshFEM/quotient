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

#include <fstream>
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
    if (supernodes[i].empty() || parents[i] == -1) {
      continue;
    }
    std::ostringstream os;
    os << "  " << parents[i] << " -> " << i << ";\n";
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

  // Extract the relevant information from the QuotientGraph.
  analysis.supernodes.resize(num_orig_vertices);
  for (Int i = 0; i < num_orig_vertices; ++i) {
    analysis.supernodes[i] = quotient_graph.FormSupernode(i);
  }
  analysis.elimination_order = quotient_graph.EliminationOrder();
  quotient_graph.ComputePostorder(&analysis.postorder);
  analysis.parents = quotient_graph.Parents();
  if (control.store_structures) {
    quotient_graph.FormEliminatedStructures(&analysis.eliminated_structures);
  }
  analysis.num_hash_collisions = quotient_graph.NumHashCollisions();
  analysis.num_hash_bucket_collisions =
      quotient_graph.NumHashBucketCollisions();
  analysis.num_aggressive_absorptions =
      quotient_graph.NumAggressiveAbsorptions();
  const std::vector<std::pair<std::string, double>>& timings =
      quotient_graph.ComponentTimes();
  for (const std::pair<std::string, double>& pairing : timings) {
    analysis.elapsed_seconds[pairing.first] = pairing.second;
  }

  return analysis;
}

} // namespace quotient

#endif // ifndef QUOTIENT_MINIMUM_DEGREE_IMPL_H_
