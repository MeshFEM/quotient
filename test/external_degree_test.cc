/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#define CATCH_CONFIG_MAIN
#include <vector>
#include "quotient.hpp"
#include "catch.hpp"

TEST_CASE("ExactExternalDegree", "[exact]") {
  using quotient::Int;
  quotient::QuotientGraph graph;
  graph.num_original_vertices = 8;
  graph.num_eliminated_vertices = 4;

  // Elimination order: 0, 7, 2.
  // Pivot: 1
  //
  // Elements: {0, 2, 7}
  // Variables: {1, 4},

  // The supernodes are:
  //   {0},
  //   {1},
  //   {2, 3},
  //   {},
  //   {4, 5, 6},
  //   {},
  //   {},
  //   {7},

  graph.supernode_sizes = std::vector<Int>{1, 1, 2, 0, 3, 0, 0, 1};

  graph.next_index = std::vector<Int>{1, 2, 3, 4, 5, 6, 7, -1};

  graph.head_index = std::vector<Int>{0, 1, 2, 2, 4, 4, 4, 7};

  graph.tail_index = std::vector<Int>{0, 1, 3, -1, 6, -1, -1, 7};

  graph.adjacency_lists = std::vector<std::vector<Int>>{
      {},
      {4},
      {},
      {},
      {1},
      {},
      {},
      {},
  };

  graph.element_lists = std::vector<std::vector<Int>>{
      {},
      {0, 2},
      {},
      {},
      {0, 7},
      {},
      {},
      {},
  };

  graph.structures = std::vector<std::vector<Int>>{
      {1, 4, 5, 6},
      {4, 5, 6},
      {1},
      {},
      {},
      {},
      {},
      {2, 3, 4, 5, 6},
  };

  graph.external_degree_heap.Reset(std::vector<Int>{
      -1, 5, -1, -1, 4, -1, -1, -1
  });
  graph.external_degree_heap.DisableIndex(0);
  graph.external_degree_heap.DisableIndex(2);
  graph.external_degree_heap.DisableIndex(3);
  graph.external_degree_heap.DisableIndex(5);
  graph.external_degree_heap.DisableIndex(6);
  graph.external_degree_heap.DisableIndex(7);

  const Int pivot = 1;
  const std::vector<Int> supernodal_pivot_structure{4};
  const std::vector<int> pivot_structure_mask{0, 0, 0, 0, 1, 0, 0, 0};
  std::vector<int> exact_degree_mask{0, 0, 0, 0, 0, 0, 0};

  const std::vector<Int> kExpectedExternalStructureSizes{
      1,  // L_0 \ L_1 = [1] = {1}
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      2,  // L_7 \ L_1 = [2] = {2, 3}
  };

  const std::vector<Int> kExpectedAggressiveAbsorptions;

  const bool aggressive_absorption = false;
  std::vector<Int> external_structure_sizes;
  std::vector<Int> aggressive_absorption_elements;
  graph.InitializeExternalStructureSizes(&external_structure_sizes);
  graph.ExternalStructureSizes(
      supernodal_pivot_structure, aggressive_absorption,
      &external_structure_sizes, &aggressive_absorption_elements);

  REQUIRE(external_structure_sizes == kExpectedExternalStructureSizes);
  REQUIRE(aggressive_absorption_elements == kExpectedAggressiveAbsorptions);

  const Int variable = 4;
  const Int exact_external_degree = ExternalDegree(
      graph, variable, pivot, pivot_structure_mask, external_structure_sizes,
      quotient::kExactExternalDegree, &exact_degree_mask);
  const Int amestoy_external_degree = ExternalDegree(
      graph, variable, pivot, pivot_structure_mask, external_structure_sizes,
      quotient::kAmestoyExternalDegree, &exact_degree_mask);
  const Int ashcraft_external_degree = ExternalDegree(
      graph, variable, pivot, pivot_structure_mask, external_structure_sizes,
      quotient::kAshcraftExternalDegree, &exact_degree_mask);
  const Int gilbert_external_degree = ExternalDegree(
      graph, variable, pivot, pivot_structure_mask, external_structure_sizes,
      quotient::kGilbertExternalDegree, &exact_degree_mask);

  // d_4 = |A_4 \ supernode(4)| + |(\cup_{e in E_4} L_e) \ supernode(4)|
  //     = |{1} \ {4, 5, 6}| + |{1, 2, 3, 4, 5, 6} \ {4, 5, 6}|
  //     = |{1}| + |{1, 2, 3}| = 1 + 3 = 4.
  REQUIRE(exact_external_degree == 4);

  // \bar{d}_4 = min(4, 4 + 0, 1 + 0 + 3) = 3.
  REQUIRE(amestoy_external_degree == 4);

  // \tilde{d}_4 = d_i if |E_i| = 2, \hat{d}_i otherwise.
  REQUIRE(ashcraft_external_degree == 4);

  // \hat{d}_4 = 1 + 3 = 4.
  REQUIRE(gilbert_external_degree == 4);
}
