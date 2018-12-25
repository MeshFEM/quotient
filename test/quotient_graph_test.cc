/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#define CATCH_CONFIG_MAIN
#include <vector>
#include "catch2/catch.hpp"

#include "quotient.hpp"

using quotient::Int;

// A reproduction of Figs. 1 and 2 from [ADD-96].
//
// The Amestoy external degree bound produces an overestimate for index
// 5 (6 with 1-based indexing) just before the third pivot is selected.
// This helps lead to a deviation between Figs. 1-2 and the AMD results.
//
TEST_CASE("ADD-96 Figures 1-2", "[ADD-96 Figs 1-2]") {
  quotient::CoordinateGraph graph;
  graph.Resize(10);

  graph.AddEdge(0, 3);
  graph.AddEdge(0, 5);
  graph.AddEdge(3, 0);
  graph.AddEdge(5, 0);

  graph.AddEdge(1, 4);
  graph.AddEdge(1, 5);
  graph.AddEdge(1, 8);
  graph.AddEdge(4, 1);
  graph.AddEdge(5, 1);
  graph.AddEdge(8, 1);

  graph.AddEdge(2, 4);
  graph.AddEdge(2, 5);
  graph.AddEdge(2, 6);
  graph.AddEdge(4, 2);
  graph.AddEdge(5, 2);
  graph.AddEdge(6, 2);

  graph.AddEdge(3, 6);
  graph.AddEdge(3, 7);
  graph.AddEdge(6, 3);
  graph.AddEdge(7, 3);

  graph.AddEdge(4, 6);
  graph.AddEdge(4, 8);
  graph.AddEdge(6, 4);
  graph.AddEdge(8, 4);

  graph.AddEdge(6, 7);
  graph.AddEdge(6, 8);
  graph.AddEdge(6, 9);
  graph.AddEdge(7, 6);
  graph.AddEdge(8, 6);
  graph.AddEdge(9, 6);

  graph.AddEdge(7, 8);
  graph.AddEdge(7, 9);
  graph.AddEdge(8, 7);
  graph.AddEdge(9, 7);

  graph.AddEdge(8, 9);
  graph.AddEdge(9, 8);

  quotient::MinimumDegreeControl control;
  control.degree_type = quotient::kExactDegree;
  control.force_minimal_pivot_indices = true;

  quotient::QuotientGraph quotient_graph(graph, control);
  const Int first_pivot = quotient_graph.FindAndProcessPivot();
  const Int first_supernode_size = quotient_graph.SupernodeSize(first_pivot);
  const std::vector<Int>& first_element = quotient_graph.Element(first_pivot);
  const Int first_cholesky_nonzeros = quotient_graph.NumPivotCholeskyNonzeros();
  const double first_cholesky_flops = quotient_graph.NumPivotCholeskyFlops();

  const Int kExpectedFirstPivot = 0;
  const Int kExpectedFirstSupernodeSize = 1;
  const std::vector<Int> kExpectedFirstElement{3, 5};
  const Int kExpectedFirstCholeskyNonzeros = 3;
  const double kExpectedFirstCholeskyFlops = 1. / 3. + 4.;

  REQUIRE(first_pivot == kExpectedFirstPivot);
  REQUIRE(first_supernode_size == kExpectedFirstSupernodeSize);
  REQUIRE(first_element == kExpectedFirstElement);
  REQUIRE(first_cholesky_nonzeros == kExpectedFirstCholeskyNonzeros);
  REQUIRE(first_cholesky_flops == kExpectedFirstCholeskyFlops);
}
