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
// Several details not discussed in the paper are worth noting:
//
//   1) The Amestoy external degree bound produces an overestimate for index
//      5 (6 with 1-based indexing) just before the third pivot is selected.
//      This helps lead to a deviation between Figs. 1-2 and the AMD results.
//
//   2) The supernode of {6, 7, 8} ({7, 8, 9} with 1-based indexing) detected
//      within G^6 in Fig. 2 of [ADD-96] is an elmination graph supernode,
//      but *not* a quotient graph supernode. Because AMD only measures the
//      latter, there is another reason for deviation.
//
// Thus, if an exact external degree computation is requested, then the
// elimination order is (0, 1, ..., 9), and the supernodes are all trivial.
//
TEST_CASE("ADD-96", "[ADD-96]") {
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

  const quotient::MinimumDegreeAnalysis analysis =
    quotient::MinimumDegree(graph, quotient::kExactExternalDegree);

  const std::vector<Int> kExpectedEliminationOrder{
      0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
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
      {6},
      {7},
      {8},
      {9},
  };

  // This structure is defined directly (modulo translation from 1-based to
  // 0-based indexing) from the bottom-right of Fig. 2 of [ADD-96].
  const std::vector<std::vector<Int>> kExpectedSupernodalStructures{
      {3, 5},
      {4, 5, 8},
      {4, 5, 6},
      {5, 6, 7},
      {5, 6, 8},
      {6, 7, 8},
      {7, 8, 9},
      {8, 9},
      {9},
      {},
  };

  REQUIRE(analysis.elimination_order == kExpectedEliminationOrder);
  REQUIRE(analysis.supernodes == kExpectedSupernodes);
  REQUIRE(analysis.supernodal_structures == kExpectedSupernodalStructures);
}
