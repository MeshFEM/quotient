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
  control.degree_type = quotient::kExactExternalDegree;
  control.force_minimal_pivot_indices = true;
  control.store_structures = true;
  control.store_aggressive_absorptions = true;
  control.store_variable_merges = true;
  const quotient::MinimumDegreeResult analysis =
      quotient::MinimumDegree(graph, control);

  const std::vector<Int> kExpectedEliminationOrder{
      0, 1, 2, 3, 4, 5, 6, 9,
  };

  // See the comment at the top of this test for why we do not expect any
  // nontrivial supernodes.
  const std::vector<std::vector<Int>> kExpectedSupernodes{
      {0},
      {1},
      {2},
      {3},
      {4},
      {5},
      {6, 7, 8},
      {},
      {},
      {9},
  };

  // This structure is defined directly (modulo translation from 1-based to
  // 0-based indexing) from the bottom-right of Fig. 2 of [ADD-96].
  const std::vector<std::vector<Int>> kExpectedPrincipalStructures{
      {3, 5},
      {4, 5, 8},
      {4, 5, 6},
      {5, 6, 7},
      {5, 6, 8},
      {6, 7, 8},
      {9},
      {},
  };

  const std::vector<std::pair<Int, Int>> kExpectedAggressiveAbsorptions;

  const std::vector<std::pair<Int, Int>> kExpectedVariableMerges{
      {6, 7},
      {6, 8},
  };

#ifdef _OPENMP
  REQUIRE(
      analysis.elimination_order.size() == kExpectedEliminationOrder.size());
  REQUIRE(
      analysis.supernodes.size() == kExpectedSupernodes.size());
  for (Int index = 0; index < 8; ++index) {
    if (index == 6) {
      // The variable merging could have picked 6, 7, or 8 as the head.
      continue;
    }
    REQUIRE(
        analysis.elimination_order[index] == kExpectedEliminationOrder[index]);
    REQUIRE(analysis.supernodes[index] == kExpectedSupernodes[index]);
  }
  REQUIRE(analysis.variable_merges.size() == kExpectedVariableMerges.size());
#else
  REQUIRE(analysis.elimination_order == kExpectedEliminationOrder);
  REQUIRE(analysis.supernodes == kExpectedSupernodes);
  REQUIRE(analysis.variable_merges == kExpectedVariableMerges);
#endif
  REQUIRE(analysis.aggressive_absorptions == kExpectedAggressiveAbsorptions);
  REQUIRE(analysis.principal_structures == kExpectedPrincipalStructures);
}

// Please see the beginning of Section 5 of [ADD-96].
TEST_CASE("ADD-96 Aggressive Absorbtion", "[ADD-96-Agg-Aborb]") {
  quotient::CoordinateGraph graph;
  graph.Resize(4);

  graph.AddEdge(0, 2);
  graph.AddEdge(0, 3);
  graph.AddEdge(2, 0);
  graph.AddEdge(3, 0);

  graph.AddEdge(1, 2);
  graph.AddEdge(1, 3);
  graph.AddEdge(2, 1);
  graph.AddEdge(3, 1);

  quotient::MinimumDegreeControl control;
  control.degree_type = quotient::kExactExternalDegree;
  control.allow_supernodes = false;
  control.force_minimal_pivot_indices = true;
  control.aggressive_absorption = true;
  control.store_structures = true;
  control.store_aggressive_absorptions = true;
  control.store_variable_merges = true;
  const quotient::MinimumDegreeResult analysis =
      quotient::MinimumDegree(graph, control);

  const std::vector<Int> kExpectedEliminationOrder{
      0, 1, 2, 3,
  };

  // See the comment at the top of this test for why we do not expect any
  // nontrivial supernodes.
  const std::vector<std::vector<Int>> kExpectedSupernodes{
      {0},
      {1},
      {2},
      {3},
  };

  // This structure is defined directly (modulo translation from 1-based to
  // 0-based indexing) from the bottom-right of Fig. 2 of [ADD-96].
  const std::vector<std::vector<Int>> kExpectedPrincipalStructures{
      {2, 3},
      {2, 3},
      {3},
      {},
  };

  // [ADD-96] discusses the aggressive absorption, 0 into 1.
  const std::vector<std::pair<Int, Int>> kExpectedAggressiveAbsorptions{
      {1, 0},
  };

  const std::vector<std::pair<Int, Int>> kExpectedVariableMerges;

  REQUIRE(analysis.elimination_order == kExpectedEliminationOrder);
  REQUIRE(analysis.supernodes == kExpectedSupernodes);
  REQUIRE(analysis.aggressive_absorptions == kExpectedAggressiveAbsorptions);
  REQUIRE(analysis.principal_structures == kExpectedPrincipalStructures);
  REQUIRE(analysis.variable_merges == kExpectedVariableMerges);
}
