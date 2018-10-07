/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>
#include "quotient.hpp"
#include "specify.hpp"

template<typename T>
double Median(const std::vector<T>& vec) {
  const std::size_t num_entries = vec.size();
  if (num_entries == 0) {
    std::cerr << "Invalid median request of empty list." << std::endl;
    return 0.;
  }

  std::vector<T> vec_copy(vec);
  std::sort(vec_copy.begin(), vec_copy.end());
  if (num_entries % 2 == 0) {
    return (vec_copy[(num_entries / 2) - 1] + vec_copy[num_entries / 2]) / 2;
  } else {
    return vec_copy[num_entries / 2];
  }
}

template<typename T>
double Mean(const std::vector<T>& vec) {
  const std::size_t num_entries = vec.size();
  if (num_entries == 0) {
    std::cerr << "Invalid mean request of empty list." << std::endl;
    return 0.;
  }

  double mean = 0.; 
  for (std::size_t index = 0; index < vec.size(); ++index) {
    mean += vec[index] / num_entries;
  }

  return mean;
}

// A helper routine for safely updating a representation of the two-norm of a
// vector as scale * sqrt(scaled_square).
void UpdateScaledSquare(double update, double* scale, double* scaled_square) {
  const double abs_update = std::abs(update);
  if (abs_update == 0.) {
    return;
  }

  if (abs_update <= *scale) {
    const double relative_scale = abs_update / *scale;
    *scaled_square += relative_scale * relative_scale;
  } else {
    const double relative_scale = *scale / abs_update;
    *scaled_square = *scaled_square * relative_scale * relative_scale + 1.;
    *scale = abs_update;
  }
}

template<typename T>
double StandardDeviation(const std::vector<T>& vec, double mean) {
  const std::size_t num_entries = vec.size();
  if (num_entries == 0) {
    std::cerr << "Invalid standard dev. request of empty list." << std::endl;
    return 0.;
  }

  double scale = 0;
  double scaled_variance = 1.;
  for (std::size_t index = 0; index < vec.size(); ++index) { 
    const double difference_from_mean = static_cast<double>(vec[index] - mean);
    UpdateScaledSquare(difference_from_mean, &scale, &scaled_variance);
  }

  return scale * std::sqrt(scaled_variance);
}

template<typename T>
void PrintMedianMeanAndStandardDeviation(
    const std::vector<T>& vec, const std::string& label) {
  const double median = Median(vec);
  const double mean = Mean(vec);
  const double sigma = StandardDeviation(vec, mean);
  std::cout << label << ": median=" << median << ", mean=" << mean 
            << ", std dev=" << sigma << std::endl;
}
 
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
  const int num_random_permutations = parser.OptionalInput<int>(
      "num_random_permutations",
      "The number of random permutations to test "
      "(in addition to the original order).",
      21);
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

  // Seed the random number generator based upon the current time.
  const unsigned srand_seed = std::time(0);
  std::cout << "Seeding std::srand with " << srand_seed << std::endl;
  std::srand(srand_seed);

  std::vector<quotient::Int> largest_supernode_sizes;
  std::vector<quotient::Int> num_strictly_lower_nonzeros;
  std::vector<double> elapsed_seconds;
  largest_supernode_sizes.reserve(num_random_permutations + 1);
  num_strictly_lower_nonzeros.reserve(num_random_permutations + 1);
  elapsed_seconds.reserve(num_random_permutations + 1); 
  for (int instance = 0; instance < num_random_permutations + 1; ++instance) {
    std::cout << "  Running analysis " << instance << " of "
              << num_random_permutations + 1 << "..." << std::endl;
    const auto start_time = std::chrono::steady_clock::now();
    quotient::MinimumDegreeAnalysis analysis = quotient::MinimumDegree(
      *graph, degree_type, allow_supernodes, aggressive_absorption,
      store_aggressive_absorptions, store_variable_merges);
    const auto end_time = std::chrono::steady_clock::now();
    const std::chrono::duration<double> duration = end_time - start_time;
    elapsed_seconds.push_back(duration.count());
    largest_supernode_sizes.push_back(analysis.LargestSupernodeSize());
    num_strictly_lower_nonzeros.push_back(analysis.NumStrictlyLowerNonzeros());
    std::cout << "  Finished analysis in " << elapsed_seconds.back()
              << " seconds. There were " << num_strictly_lower_nonzeros.back()
              << " subdiagonal nonzeros and the largest supernode had "
              << largest_supernode_sizes.back() << " members." << std::endl;
    if (instance == num_random_permutations) {
      break;
    }

    // Generate a random permutation.
    std::vector<quotient::Int> permutation(graph->NumSources());
    std::iota(permutation.begin(), permutation.end(), 0);
    std::random_shuffle(permutation.begin(), permutation.end());

    // Apply the permutation to the graph.
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

  PrintMedianMeanAndStandardDeviation(
      num_strictly_lower_nonzeros, "Num strictly lower nonzeros");
  PrintMedianMeanAndStandardDeviation(
      largest_supernode_sizes, "Largest supernode sizes");
  PrintMedianMeanAndStandardDeviation(elapsed_seconds, "Elapsed seconds");

  return 0;
}
