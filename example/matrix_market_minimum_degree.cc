/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <numeric>
#include <vector>
#include "quotient.hpp"
#include "specify.hpp"

// A simple data structure for storing the median, mean, and standard deviation
// of a random variable.
struct SufficientStatistics {
  // The median of the set of samples.
  double median;

  // The mean of the set of samples.
  double mean;

  // The standard deviation of the set of samples.
  double standard_deviation;
};

// Pretty prints the SufficientStatistics structure.
void PrintSufficientStatistics(
    const SufficientStatistics& stats, const std::string& label) {
  std::cout << label << ": median=" << stats.median << ", mean=" << stats.mean
            << ", stddev=" << stats.standard_deviation << std::endl;
}

// A list of properties to measure from an AMD reordering.
struct AMDExperiment {
  // The number of (structural) nonzeros in the strictly lower triangle of the
  // associated Cholesky factor.
  SufficientStatistics num_lower_nonzeros;

  // The number of members of the largest supernode.
  SufficientStatistics largest_supernode_size;

  // The number of seconds that elapsed during the AMD analysis.
  SufficientStatistics elapsed_seconds;
};

// Pretty prints the AMDExperiment structure.
void PrintAMDExperiment(
    const AMDExperiment& experiment, const std::string& label) {
  std::cout << label << ":\n";
  PrintSufficientStatistics(
      experiment.num_lower_nonzeros, "  num_lower_nonzeros");
  PrintSufficientStatistics(
      experiment.largest_supernode_size, "  largest_supernode_size");
  PrintSufficientStatistics(
      experiment.elapsed_seconds, "  elapsed_seconds");
}

// Returns the median of a given vector.
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

// Returns the mean of a given vector.
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

// Returns the standard deviation of a given vector with a given mean.
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

// Returns the SufficientStatistics properties for a given vector.
template<typename T>
SufficientStatistics GetSufficientStatistics(const std::vector<T>& vec) {
  SufficientStatistics stats;
  stats.median = Median(vec);
  stats.mean = Mean(vec);
  stats.standard_deviation = StandardDeviation(vec, stats.mean);
  return stats;
}

// Returns the AMDExperiment statistics for a single Matrix Market input matrix.
AMDExperiment RunMatrixMarketAMDTest(
    const std::string& filename,
    const quotient::MinimumDegreeControl& control,
    bool force_symmetry,
    int num_random_permutations,
    bool print_progress) {
  if (print_progress) {
    std::cout << "Reading CoordinateGraph from " << filename << "..."
              << std::endl;
  }
  std::unique_ptr<quotient::CoordinateGraph> graph =
      quotient::CoordinateGraph::FromMatrixMarket(filename);
  if (!graph) {
    std::cerr << "Could not open " << filename << "." << std::endl;
    AMDExperiment experiment;
    return experiment;
  }
  if (print_progress) {
    std::cout << "Graph had " << graph->NumSources() << " sources and "
              << graph->NumEdges() << " edges." << std::endl;
  }

  // Force symmetry since many of the examples are not. We form the nonzero
  // pattern of A + A'.
  if (force_symmetry) {
    if (print_progress) {
      std::cout << "Enforcing graph symmetry..." << std::endl;
    }
    graph->ReserveEdgeAdditions(graph->NumEdges());
    for (const std::pair<quotient::Int, quotient::Int>& edge : graph->Edges()) {
      graph->QueueEdgeAddition(edge.second, edge.first);
    }
    graph->FlushEdgeQueues();
  }

  std::vector<quotient::Int> largest_supernode_sizes;
  std::vector<quotient::Int> num_strictly_lower_nonzeros;
  std::vector<double> elapsed_seconds;
  largest_supernode_sizes.reserve(num_random_permutations + 1);
  num_strictly_lower_nonzeros.reserve(num_random_permutations + 1);
  elapsed_seconds.reserve(num_random_permutations + 1); 
  for (int instance = 0; instance < num_random_permutations + 1; ++instance) {
    if (print_progress) {
      std::cout << "  Running analysis " << instance << " of "
                << num_random_permutations + 1 << "..." << std::endl;
    }
    quotient::Timer timer;
    timer.Start();
    const quotient::MinimumDegreeAnalysis analysis = quotient::MinimumDegree(
      *graph, control);
    elapsed_seconds.push_back(timer.Stop());
    largest_supernode_sizes.push_back(analysis.LargestSupernodeSize());
    num_strictly_lower_nonzeros.push_back(analysis.NumStrictlyLowerNonzeros());
    if (print_progress) {
      std::cout << "  Finished analysis in " << elapsed_seconds.back()
                << " seconds. There were " << num_strictly_lower_nonzeros.back()
                << " subdiagonal nonzeros and the largest supernode had "
                << largest_supernode_sizes.back() << " members." << std::endl;
    }
    if (control.store_pivot_element_list_sizes) {
      std::cout << "  Fraction of pivots with multiple elements: "
                << analysis.FractionOfPivotsWithMultipleElements() << std::endl;
    }
    if (control.store_num_degree_updates_with_multiple_elements) {
      std::cout << "  Fraction of degree updates with multiple elements: "
                << analysis.FractionOfDegreeUpdatesWithMultipleElements()
                << std::endl;
    }
    if (control.time_stages) {
      for (const std::pair<std::string, double>& pairing :
          analysis.elapsed_seconds) {
        std::cout << "    " << pairing.first << ": " << pairing.second
                  << " seconds." << std::endl;
      }
    }
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

  AMDExperiment experiment;
  experiment.num_lower_nonzeros = GetSufficientStatistics(
      num_strictly_lower_nonzeros);
  experiment.largest_supernode_size = GetSufficientStatistics(
      largest_supernode_sizes);
  experiment.elapsed_seconds = GetSufficientStatistics(elapsed_seconds);

  return experiment;
}

// Returns a map from the identifying string of each test matrix from the
// Amestoy/Davis/Duff Approximate Minimum Degree reordering 1996 paper meant
// to loosely reproduce Fig. 2.
//
// It is worth noting that the LHR34 results from the paper appear to be
// incorrect, as the results shown in 
//
//   https://www.cise.ufl.edu/research/sparse/matrices/Mallya/lhr34.html
//
// agree with the results observed from this code's implementation.
//
std::unordered_map<std::string, AMDExperiment> RunADD96Tests(
    const std::string& matrix_market_directory,
    const quotient::MinimumDegreeControl& control,
    int num_random_permutations,
    bool print_progress) {
  const std::vector<std::string> kMatrixNames{
      "appu",
      "bbmat",
      "bcsstk30",
      "bcsstk31",
      "bcsstk32",
      "bcsstk33",
      "crystk02",
      "crystk03",
      "ct20stif",
      "ex11",
      "ex19",
      "ex40",
      "finan512",
      "lhr34",
      "lhr71",
      "nasasrb",
      "olafu",
      "orani678",
      "psmigr_1",
      "raefsky1",
      "raefsky3",
      "raefsky4",
      "rim",
      "venkat01",
      "wang3",
      "wang4",
  };
  const bool force_symmetry = true;

  std::unordered_map<std::string, AMDExperiment> experiments;
  for (const std::string& matrix_name : kMatrixNames) {
    const std::string filename = matrix_market_directory + "/" + matrix_name +
        "/" + matrix_name + ".mtx";
    experiments[matrix_name] = RunMatrixMarketAMDTest(
        filename,
        control,
        force_symmetry,
        num_random_permutations,
        print_progress);
  }

  return experiments;
}

int main(int argc, char** argv) {
  specify::ArgumentParser parser(argc, argv);
  const std::string filename = parser.OptionalInput<std::string>(
      "filename", "The location of a Matrix Market file.", "");
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
      false);
  const bool store_variable_merges = parser.OptionalInput<bool>(
      "store_variable_merges",
      "Store the variable merge list?",
      false);
  const bool store_pivot_element_list_sizes =
      parser.OptionalInput<bool>(
          "store_pivot_element_list_sizes",
          "Store the length of each pivot's element list?",
          false);
  const bool store_num_degree_updates_with_multiple_elements =
      parser.OptionalInput<bool>(
          "store_num_degree_updates_with_multiple_elements",
          "Store the number of degree updates whose corresponding variable had "
          "more than two members in its element list?",
          false);
  const int num_random_permutations = parser.OptionalInput<int>(
      "num_random_permutations",
      "The number of random permutations to test "
      "(in addition to the original order).",
      21);
  const bool force_symmetry = parser.OptionalInput<bool>(
      "force_symmetry",
      "Use the nonzero pattern of A + A'?",
      true);
  const bool print_progress = parser.OptionalInput<bool>(
      "print_progress",
      "Print the progress of the experiments?",
      false);
  const std::string matrix_market_directory = parser.OptionalInput<std::string>(
      "matrix_market_directory",
      "The directory where the ADD96 matrix market .tar.gz's were unpacked",
      "");
  const bool time_stages = parser.OptionalInput<bool>(
      "time_stages",
      "Report the timings of each stage of MinimumDegree?",
      false);
  if (!parser.OK()) {
    return 0;
  }
  if (filename.empty() && matrix_market_directory.empty()) {
    std::cerr << "One of 'filename' or 'matrix_market_directory' must be "
                 "specified.\n" << std::endl;
    parser.PrintReport();
    return 0;
  }

  quotient::MinimumDegreeControl control;
  control.degree_type =
      static_cast<quotient::ExternalDegreeType>(degree_type_int);
  control.allow_supernodes = allow_supernodes;
  control.aggressive_absorption = aggressive_absorption;
  control.store_aggressive_absorptions = store_aggressive_absorptions;
  control.store_variable_merges = store_variable_merges;
  control.store_pivot_element_list_sizes = store_pivot_element_list_sizes;
  control.store_num_degree_updates_with_multiple_elements =
      store_num_degree_updates_with_multiple_elements;
  control.time_stages = time_stages;

  // Seed the random number generator based upon the current time.
  const unsigned srand_seed = std::time(0);
  std::cout << "Seeding std::srand with " << srand_seed << std::endl;
  std::srand(srand_seed);

  if (!matrix_market_directory.empty()) {
    const std::unordered_map<std::string, AMDExperiment> experiments =
        RunADD96Tests(
            matrix_market_directory,
            control,
            num_random_permutations,
            print_progress);
    for (const std::pair<std::string, AMDExperiment>& pairing : experiments) {
      PrintAMDExperiment(pairing.second, pairing.first);  
    }
  } else {
    const AMDExperiment experiment = RunMatrixMarketAMDTest(
        filename,
        control,
        force_symmetry,
        num_random_permutations,
        print_progress);
    PrintAMDExperiment(experiment, filename);
  }

  return 0;
}
