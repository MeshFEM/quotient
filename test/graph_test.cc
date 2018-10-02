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

TEST_CASE("Simple example", "[simple]") {
  quotient::CoordinateGraph graph;
  graph.Resize(5);
  graph.AddEdge(0, 1);
  graph.AddEdge(3, 2);
  graph.AddEdge(0, 4);
  graph.AddEdge(2, 0);
  graph.AddEdge(3, 2);
  graph.AddEdge(0, 2);

  // The sorted edges should now be:
  //   (0, 1), (0, 2), (0, 4), (2, 0), (3, 2).

  const std::vector<quotient::GraphEdge>& edges = graph.Edges();
  REQUIRE(edges.size() == 5);
  REQUIRE(edges[0] == quotient::GraphEdge(0, 1));
  REQUIRE(edges[1] == quotient::GraphEdge(0, 2));
  REQUIRE(edges[2] == quotient::GraphEdge(0, 4));
  REQUIRE(edges[3] == quotient::GraphEdge(2, 0));
  REQUIRE(edges[4] == quotient::GraphEdge(3, 2));

  graph.RemoveEdge(2, 0);
  graph.RemoveEdge(0, 2);
  graph.RemoveEdge(0, 2);
  
  REQUIRE(edges.size() == 3);
  REQUIRE(edges[0] == quotient::GraphEdge(0, 1));
  REQUIRE(edges[1] == quotient::GraphEdge(0, 4));
  REQUIRE(edges[2] == quotient::GraphEdge(3, 2));

  graph.AddEdge(0, 0);

  REQUIRE(edges.size() == 4);
  REQUIRE(edges[0] == quotient::GraphEdge(0, 0));
  REQUIRE(edges[1] == quotient::GraphEdge(0, 1));
  REQUIRE(edges[2] == quotient::GraphEdge(0, 4));
  REQUIRE(edges[3] == quotient::GraphEdge(3, 2));
  REQUIRE(graph.EdgeExists(0, 4));
  REQUIRE(!graph.EdgeExists(0, 3));
  REQUIRE(graph.NumConnections(0) == 3);
  REQUIRE(graph.NumConnections(3) == 1);
  REQUIRE(graph.SourceEdgeOffset(0) == 0);
  REQUIRE(graph.SourceEdgeOffset(1) == 3);
  REQUIRE(graph.SourceEdgeOffset(2) == 3);
  REQUIRE(graph.SourceEdgeOffset(3) == 3);
  REQUIRE(graph.EdgeOffset(1, 1) == 3);
  REQUIRE(graph.EdgeOffset(3, 2) == 3);
  REQUIRE(graph.EdgeOffset(3, 3) == 4);
}

TEST_CASE("Batch example", "[batch]") {
  quotient::CoordinateGraph graph;
  graph.Resize(5);

  graph.ReserveEdgeAdditions(6);
  graph.QueueEdgeAddition(0, 1);
  graph.QueueEdgeAddition(3, 2);
  graph.QueueEdgeAddition(0, 4);
  graph.QueueEdgeAddition(2, 0);
  graph.QueueEdgeAddition(3, 2);
  graph.QueueEdgeAddition(0, 2);
  graph.FlushEdgeQueues();

  graph.ReserveEdgeRemovals(3);
  graph.QueueEdgeRemoval(2, 0);
  graph.QueueEdgeRemoval(0, 2);
  graph.QueueEdgeRemoval(0, 2);
  graph.FlushEdgeQueues();

  const std::vector<quotient::GraphEdge>& edges = graph.Edges();
  
  REQUIRE(edges.size() == 3);
  REQUIRE(edges[0] == quotient::GraphEdge(0, 1));
  REQUIRE(edges[1] == quotient::GraphEdge(0, 4));
  REQUIRE(edges[2] == quotient::GraphEdge(3, 2));

  graph.ReserveEdgeAdditions(2);
  graph.QueueEdgeAddition(3, 3);
  graph.QueueEdgeAddition(2, 2);
  graph.FlushEdgeQueues();

  REQUIRE(edges.size() == 5);
  REQUIRE(edges[0] == quotient::GraphEdge(0, 1));
  REQUIRE(edges[1] == quotient::GraphEdge(0, 4));
  REQUIRE(edges[2] == quotient::GraphEdge(2, 2));
  REQUIRE(edges[3] == quotient::GraphEdge(3, 2));
  REQUIRE(edges[4] == quotient::GraphEdge(3, 3));
}
