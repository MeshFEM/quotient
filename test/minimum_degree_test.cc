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
  control.degree_type = quotient::kExactDegree;
  control.force_minimal_pivot_indices = true;
  control.store_supernodes = true;
  control.store_structures = true;
  quotient::MinimumDegreeResult analysis =
      quotient::MinimumDegree(graph, control);
  for (std::size_t i = 0; i < analysis.supernodes.size(); ++i) {
    std::sort(analysis.supernodes[i].begin(), analysis.supernodes[i].end());
  }
  for (std::size_t i = 0; i < analysis.eliminated_structures.size(); ++i) {
    std::sort(
        analysis.eliminated_structures[i].begin(),
        analysis.eliminated_structures[i].end());
  }

  // Because of the ordering of the hash bucket, we will prefer the last member
  // as the key. But there are several equally-valid solutions (the
  // only nontrivial supervariable should be {6, 7, 8}, but the principal
  // member can vary).
  const std::vector<Int> kExpectedEliminationOrder{
      0, 1, 2, 3, 4, 5, 8, 9,
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
      {},
      {},
      {6, 7, 8},
      {9},
  };

  // This structure is defined directly (modulo translation from 1-based to
  // 0-based indexing) from the bottom-right of Fig. 2 of [ADD-96].
  const std::vector<std::vector<Int>> kExpectedEliminatedStructures{
      {3, 5},
      {4, 5, 8},
      {4, 5, 6},
      {5, 6, 7},
      {5, 6, 8},
      {6, 7, 8},
      {9},
      {},
  };

  const Int kExpectedNumAggressiveAbsorptions = 0;

  REQUIRE(analysis.elimination_order == kExpectedEliminationOrder);
  REQUIRE(analysis.supernodes == kExpectedSupernodes);
  REQUIRE(
      analysis.num_aggressive_absorptions == kExpectedNumAggressiveAbsorptions);
  REQUIRE(analysis.eliminated_structures == kExpectedEliminatedStructures);
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
  control.degree_type = quotient::kExactDegree;
  control.allow_supernodes = false;
  control.force_minimal_pivot_indices = true;
  control.aggressive_absorption = true;
  control.store_supernodes = true;
  control.store_structures = true;
  quotient::MinimumDegreeResult analysis =
      quotient::MinimumDegree(graph, control);
  for (std::size_t i = 0; i < analysis.supernodes.size(); ++i) {
    std::sort(analysis.supernodes[i].begin(), analysis.supernodes[i].end());
  }
  for (std::size_t i = 0; i < analysis.eliminated_structures.size(); ++i) {
    std::sort(
        analysis.eliminated_structures[i].begin(),
        analysis.eliminated_structures[i].end());
  }

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
  const std::vector<std::vector<Int>> kExpectedEliminatedStructures{
      {2, 3},
      {2, 3},
      {3},
      {},
  };

  // [ADD-96] discusses the aggressive absorption, 0 into 1.
  const Int kExpectedNumAggressiveAbsorptions = 1;

  REQUIRE(analysis.elimination_order == kExpectedEliminationOrder);
  REQUIRE(analysis.supernodes == kExpectedSupernodes);
  REQUIRE(
      analysis.num_aggressive_absorptions == kExpectedNumAggressiveAbsorptions);
  REQUIRE(analysis.eliminated_structures == kExpectedEliminatedStructures);
}
