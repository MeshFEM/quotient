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

TEST_CASE("VariablesAreQuotientIndistinguishable", "[indistinguishable]") {
  quotient::QuotientGraph graph;
  graph.adjacency_lists = std::vector<std::vector<quotient::Int>>{
      {1, 5, 7, 15, 27},
      {0, 5, 7, 15, 27},
  };

  graph.element_lists = std::vector<std::vector<quotient::Int>>{
      {},
      {0},
  };

  REQUIRE(graph.VariablesAreQuotientIndistinguishable(0, 1) == true);
}
