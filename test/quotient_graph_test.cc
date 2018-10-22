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

TEST_CASE("StructuralSupervariablesAreQuotientIndistinguishable",
          "[indistinguishable]") {
  quotient::QuotientGraph graph;
  graph.adjacency_lists = std::vector<std::vector<quotient::Int>>{
      {5, 7, 15, 27},
      {5, 7, 15, 27},
  };

  graph.element_lists = std::vector<std::vector<quotient::Int>>{
      {0},
      {0},
  };

  REQUIRE(graph.StructuralSupervariablesAreQuotientIndistinguishable(0, 0) ==
      true);
  REQUIRE(graph.StructuralSupervariablesAreQuotientIndistinguishable(0, 1) ==
      true);
  REQUIRE(graph.StructuralSupervariablesAreQuotientIndistinguishable(1, 0) ==
      true);
  REQUIRE(graph.StructuralSupervariablesAreQuotientIndistinguishable(1, 1) ==
      true);

  graph.adjacency_lists = std::vector<std::vector<quotient::Int>>{
      {5, 8, 15, 27},
      {5, 7, 15, 27},
  };

  REQUIRE(graph.StructuralSupervariablesAreQuotientIndistinguishable(0, 0) ==
      true);
  REQUIRE(graph.StructuralSupervariablesAreQuotientIndistinguishable(0, 1) ==
      false);
  REQUIRE(graph.StructuralSupervariablesAreQuotientIndistinguishable(1, 0) ==
      false);
  REQUIRE(graph.StructuralSupervariablesAreQuotientIndistinguishable(1, 1) ==
      true);

  graph.adjacency_lists = std::vector<std::vector<quotient::Int>>{
      {},
      {30301, 30535},
  };

  graph.element_lists = std::vector<std::vector<quotient::Int>>{
      {29827, 29965, 30304, 30532},
      {29827, 29965},
  };

  REQUIRE(graph.StructuralSupervariablesAreQuotientIndistinguishable(0, 0) ==
      true);
  REQUIRE(graph.StructuralSupervariablesAreQuotientIndistinguishable(0, 1) ==
      false);
  REQUIRE(graph.StructuralSupervariablesAreQuotientIndistinguishable(1, 0) ==
      false);
  REQUIRE(graph.StructuralSupervariablesAreQuotientIndistinguishable(1, 1) ==
      true);
}
