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
  quotient::QuotientGraph graph;
  graph.num_original_vertices = 8;
  graph.num_eliminated_vertices = 4;

  // Elimination order: 0, 7, 2.
  // Pivot: 1
  //
  // Elements: {0, 2, 7}
  // Variables: {1, 4},

  graph.supernodes = std::vector<std::vector<quotient::Int>>{
      {0},
      {1},
      {2, 3},
      {},
      {4, 5, 6},
      {},
      {},
      {7},
  };

  graph.adjacency_lists = std::vector<std::vector<quotient::Int>>{
      {},
      {4, 6},
      {},
      {},
      {1},
      {},
      {},
      {},
  };

  graph.element_lists = std::vector<std::vector<quotient::Int>>{
      {},
      {0, 2},
      {},
      {},
      {0, 7},
      {},
      {},
      {},
  };

  graph.structures = std::vector<std::vector<quotient::Int>>{
      {1, 4, 5, 6},
      {4, 5, 6},
      {1},
      {},
      {},
      {},
      {},
      {3, 4, 5, 6},
  };

  graph.external_degree_heap.Reset(std::vector<quotient::Int>{
      -1, 5, -1, -1, 3, -1, -1, -1
  });
  graph.external_degree_heap.DisableIndex(0);
  graph.external_degree_heap.DisableIndex(2);
  graph.external_degree_heap.DisableIndex(3);
  graph.external_degree_heap.DisableIndex(5);
  graph.external_degree_heap.DisableIndex(6);
  graph.external_degree_heap.DisableIndex(7);

  const quotient::Int pivot = 1;

  // '4' is the only principal member of L_1 = {4, 5, 6}.
  const std::vector<quotient::Int> kExpectedSupernodalPivotStructure{4};

  const std::vector<quotient::Int> supernodal_pivot_structure =
      graph.FormSupernodalStructure(pivot);

  REQUIRE(supernodal_pivot_structure == kExpectedSupernodalPivotStructure);

  const std::vector<quotient::Int> kExpectedExternalStructureSizes{
      1,  // L_0 \ L_1 = {1}
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      1,  // L_7 \ L_1 = {3}
  };

  const std::vector<quotient::Int> kExpectedAggressiveAbsorptions;

  const bool aggressive_absorption = false;
  std::vector<quotient::Int> external_structure_sizes;
  std::vector<quotient::Int> aggressive_absorption_elements;
  graph.InitializeExternalStructureSizes(&external_structure_sizes);
  graph.ExternalStructureSizes(
      supernodal_pivot_structure, aggressive_absorption,
      &external_structure_sizes, &aggressive_absorption_elements);

  REQUIRE(external_structure_sizes == kExpectedExternalStructureSizes);
  REQUIRE(aggressive_absorption_elements == kExpectedAggressiveAbsorptions);

  const quotient::Int variable = 4;
  const quotient::Int exact_external_degree = ExternalDegree(
      graph, variable, pivot, external_structure_sizes,
      quotient::kExactExternalDegree);
  const quotient::Int amestoy_external_degree = ExternalDegree(
      graph, variable, pivot, external_structure_sizes,
      quotient::kAmestoyExternalDegree);
  const quotient::Int ashcraft_external_degree = ExternalDegree(
      graph, variable, pivot, external_structure_sizes,
      quotient::kAshcraftExternalDegree);
  const quotient::Int gilbert_external_degree = ExternalDegree(
      graph, variable, pivot, external_structure_sizes,
      quotient::kGilbertExternalDegree);

  // d_4 = |A_4 \ supernode(4)| + |(\cup_{e in E_4} L_e) \ supernode(4)|
  //     = |{1} \ {4, 5, 6}| + ||({1, 4, 5} \cup {3, 5}) \ {4, 5, 6}|
  //     = 1 + 2 = 3.
  REQUIRE(exact_external_degree == 3);

  // \bar{d}_4 = min(4, 3 + 0, 1 + 0 + 2) = 3.
  REQUIRE(amestoy_external_degree == 3);

  // \tilde{d}_4 = d_i if |E_i| = 2, \hat{d}_i otherwise.
  REQUIRE(ashcraft_external_degree == 3);

  // \hat{d}_4 = 1 + 2 = 3.
  REQUIRE(gilbert_external_degree == 3);
}
