/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>
#include "quotient.hpp"
#include "specify.hpp"
 
int main(int argc, char** argv) {
  specify::ArgumentParser parser(argc, argv);
  const std::string filename = parser.RequiredInput<std::string>(
      "filename", "The location of a Matrix Market file.");
  const int degree_type_int = parser.OptionalInput<int>(
      "degree_type_int",
      "The degree approximation type.\n"
      "0:exact, 1:Amestoy, 2:Ashcraft, 3:Gilbert",
      1);
  const bool allow_supernodes = parser.OptionalInput<bool>(
      "allow_supernodes",
      "Allow variables to be merged into supernodes?",
      true);
  const bool aggressive_absorption = parser.OptionalInput<bool>(
      "aggressive_absorption",
      "Eliminate elements with aggressive absorption?",
      true);
  const bool store_aggressive_absorptions = parser.OptionalInput<bool>(
      "store_aggressive_absorptions",
      "Store the aggressive absorption list?",
      true);
  const bool store_variable_merges = parser.OptionalInput<bool>(
      "store_variable_merges",
      "Store the variable merge list?",
      true);
  const bool randomly_permute = parser.OptionalInput<bool>(
      "randomly_permute",
      "randomly permute initial ordering?",
      false);
  const bool force_symmetry = parser.OptionalInput<bool>(
      "force_symmetry",
      "Use the nonzero pattern of A + A'?",
      true);
  if (!parser.OK()) {
    return 0;
  }

  const quotient::ExternalDegreeType degree_type =
      static_cast<quotient::ExternalDegreeType>(degree_type_int);

  std::cout << "Reading CoordinateGraph from " << filename << "..."
            << std::endl;
  std::unique_ptr<quotient::CoordinateGraph> graph =
      quotient::CoordinateGraph::FromMatrixMarket(filename);
  if (!graph) {
    std::cerr << "Could not open " << filename << "." << std::endl;
    return 0;
  }
  std::cout << "Graph had " << graph->NumSources() << " sources and "
            << graph->NumEdges() << " edges." << std::endl;

  if (randomly_permute) {
    // Seed the random number generator based upon the current time.
    std::srand(unsigned(std::time(0)));

    std::vector<quotient::Int> permutation(graph->NumSources());
    std::iota(permutation.begin(), permutation.end(), 0);
    std::random_shuffle(permutation.begin(), permutation.end());

    quotient::CoordinateGraph permuted_graph;
    permuted_graph.Resize(graph->NumSources());
    permuted_graph.ReserveEdgeAdditions(graph->NumEdges());
    for (const std::pair<quotient::Int, quotient::Int>& edge : graph->Edges()) {
      permuted_graph.QueueEdgeAddition(
          permutation[edge.first], permutation[edge.second]);
    }
    permuted_graph.FlushEdgeQueues();

    *graph = permuted_graph;
  }

  // Force symmetry since many of the examples are not. We form the nonzero
  // pattern of A + A'.
  if (force_symmetry) {
    std::cout << "Enforcing graph symmetry..." << std::endl;
    graph->ReserveEdgeAdditions(graph->NumEdges());
    for (const std::pair<quotient::Int, quotient::Int>& edge : graph->Edges()) {
      graph->QueueEdgeAddition(edge.second, edge.first);
    }
    graph->FlushEdgeQueues();
  }

  std::cout << "Running MinimumDegree analysis..." << std::endl;
  quotient::MinimumDegreeAnalysis analysis = quotient::MinimumDegree(
    *graph, degree_type, allow_supernodes, aggressive_absorption,
    store_aggressive_absorptions, store_variable_merges);
  std::cout << "Finished MinimumDegree." << std::endl;

  const quotient::Int largest_supernode = analysis.LargestSupernode();
  const quotient::Int largest_supernode_size =
      analysis.supernodes[largest_supernode].size();
  std::cout << "Largest supernode: " << largest_supernode << ", "
            << largest_supernode_size << " entries." << std::endl;

  std::cout << "Num strictly-lower nonzeros: "
            << analysis.NumStrictlyLowerNonzeros() << std::endl;

  return 0;
}
