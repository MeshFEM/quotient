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
#include "quotient/timer.hpp"

#include "quotient/minimum_degree.hpp"

namespace quotient {

inline void ForestToDot(const std::string& filename,
                        const Buffer<Int>& parents) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Could not open " << filename << std::endl;
    return;
  }

  file << "digraph g{\n";
  for (std::size_t i = 0; i < parents.Size(); ++i) {
    if (parents[i] < 0) {
      continue;
    }
    std::ostringstream os;
    os << "  " << parents[i] << " -> " << i << ";\n";
    file << os.str();
  }
  file << "}\n";
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

inline MinimumDegreeResult MinimumDegree(QuotientGraph* graph) {
  const Int num_vertices = graph->NumVertices();
  const MinimumDegreeControl& control = graph->Control();

  // Eliminate the variables.
  MinimumDegreeResult analysis;
  if (control.store_pivot_element_list_sizes) {
    analysis.pivot_element_list_sizes.reserve(num_vertices);
  }
  if (control.store_num_degree_updates_with_multiple_elements) {
    analysis.num_degree_updates_with_multiple_elements = 0;
  }
  while (graph->NumEliminatedVertices() < num_vertices) {
    graph->FindAndProcessPivot();

    if (control.store_pivot_element_list_sizes) {
      analysis.pivot_element_list_sizes.push_back(graph->NumPivotElements());
    }
    analysis.num_cholesky_nonzeros += graph->NumPivotCholeskyNonzeros();
    analysis.num_cholesky_flops += graph->NumPivotCholeskyFlops();
    analysis.num_degree_updates += graph->NumPivotDegreeUpdates();
    if (control.store_num_degree_updates_with_multiple_elements) {
      analysis.num_degree_updates_with_multiple_elements +=
          graph->NumPivotDegreeUpdatesWithMultipleElements();
    }
  }
  graph->CombineDenseNodes();

  // Assume the Schur complement of the non-dense supernodes onto the "dense"
  // ones results in a dense Schur complement.
  const Int num_dense = graph->NumDense();
  analysis.num_cholesky_nonzeros += ((num_dense + 1) * num_dense) / 2;
  analysis.num_cholesky_flops += std::pow(1. * num_dense, 3.) / 3.;

  analysis.num_hash_collisions = graph->NumHashCollisions();
  analysis.num_hash_bucket_collisions = graph->NumHashBucketCollisions();
  analysis.num_aggressive_absorptions = graph->NumAggressiveAbsorptions();

  return analysis;
}

}  // namespace quotient

#endif  // ifndef QUOTIENT_MINIMUM_DEGREE_IMPL_H_
