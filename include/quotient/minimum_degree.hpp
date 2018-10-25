/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_MINIMUM_DEGREE_H_
#define QUOTIENT_MINIMUM_DEGREE_H_

#include <unordered_map>
#include <vector>

#include "quotient/config.hpp"
#include "quotient/coordinate_graph.hpp"
#include "quotient/external_degree.hpp"

namespace quotient {

// A data structure for controlling the MinimumDegree reordering routine.
struct MinimumDegreeControl {
  // The type of approximation to use for the external degree estimates.
  ExternalDegreeType degree_type = kAmestoyExternalDegree;

  // The type of hash function to use for converting a variable to std::size_t.
  VariableHashType hash_type = kBasicVariableHash;

  // Whether nontrivial supernodes are allowed. It is highly recommended that
  // this remain true.
  bool allow_supernodes = true;

  // Whether aggressive element absorptions are allowed.
  bool aggressive_absorption = true;

  // Whether the entire degree list will be traversed in order to ensure that
  // the member of the minimal degree list with smallest index is chosen.
  bool force_minimal_pivot_indices = false;

  // Whether the list of pairs of aggressive element absorptions should be
  // returned in the MinimumDegreeAnalysis result of MinimumDegree.
  bool store_aggressive_absorptions = false;

  // Whether the list of pairs of variable merges should be returned in the
  // MinimumDegreeAnalysis result of MinimumDegree.
  bool store_variable_merges = false;

  // Whether a list should be stored of the lengths of the element lists of
  // the pivots.
  bool store_pivot_element_list_sizes = false;

  // Whether or not to store the count of the number of external degree updates
  // which involved more than two (and, separately, how many less than or equal
  // to two) elements in the variable's element list.
  bool store_num_degree_updates_with_multiple_elements = false;

  // Whether or not to create and store the nonzeros structures of each pivot.
  bool store_structures = false;

  // Whether a breakdown of the elapsed seconds of each stage of the reordering
  // should be saved.
  bool time_stages = false;
};


// The result of running the MinimumDegree reordering algorithm. It contains
// the ordered list of eliminated principal vertices, the list of supernodes,
// and the supernodal nonzero structure of each principal column.
struct MinimumDegreeAnalysis {
  // The recommended elimination order of the supernodes (with each supernode
  // represented by its principal member).
  std::vector<Int> elimination_order;

  // The list of supernodes. Entry 'i' is nonempty if and only if 'i' is the
  // principal member of the supernode, in which case 'supernodes[i]' contains
  // the sorted list of members of the supernode.
  std::vector<std::vector<Int>> supernodes;

  // The structures of the supernodes. Entry 'index' corresponds to the
  // structure of supernode 'elimination_order[index]'.
  std::vector<std::vector<Int>> principal_structures;

  // An optional list (based on the value of
  // 'MinimumDegreeControl.store_aggressive_absorptions') of aggressive element
  // absorption pairs: each pair (e, f) consists of the absorbing element, e,
  // and the absorbed element, f.
  std::vector<std::pair<Int, Int>> aggressive_absorptions;

  // An optional list (based on the value of
  // 'MinimumDegreeControl.store_variable_merges') of supervariable merge
  // pairs: each pair (i, j) consists of the absorbing supervariable, i, and
  // the absorbed supervariable, j.
  std::vector<std::pair<Int, Int>> variable_merges;

  // An optional list (based on the value of
  // 'MinimumDegreeControl.store_pivot_element_list_sizes') of the lengths of
  // the element lists of each pivot.
  std::vector<Int> pivot_element_list_sizes;

  // An optional count (based on the value of
  // 'MinimumDegreeControl.store_num_degree_updates_with_multiple_elements')
  // of the number of external degree updates that involved a variable with
  // more than two members in its element list.
  Int num_degree_updates_with_multiple_elements = -1;

  // An optional count (based on the value of
  // 'MinimumDegreeControl.store_num_degree_updates_with_multiple_elements')
  // of the number of external degree updates that involved a variable with
  // more than two members in its element list.
  Int num_degree_updates_without_multiple_elements = -1;

  // We will push to elimination order as the reordering algorithm progresses,
  // so we will allocate an upper bound for the amount of required space.
  // The 'supernodes' and 'structures' variables will be copied over from
  // the quotient graph just before the analysis completes.
  MinimumDegreeAnalysis(Int num_vertices);

  // Returns the number of structural nonzeros in the strictly lower-triangular
  // factor.
  Int NumStrictlyLowerNonzeros() const;

  // Returns the principal member of the largest supernode.
  Int LargestSupernode() const;

  // Returns the number of members of the largest supernode.
  Int LargestSupernodeSize() const;

  // Returns the fraction of pivots whose element list had more than two
  // members.
  double FractionOfPivotsWithMultipleElements() const;

  // Returns the fraction of degree updates whose corresponding variable had an
  // element list had more than two members.
  double FractionOfDegreeUpdatesWithMultipleElements() const;

  // An optional (based upon the value of 'MinimumDegreeControl.time_stages')
  // map from the stage names to the corresponding elapsed seconds.
  std::unordered_map<std::string, double> elapsed_seconds;
};

// Returns a supernodal reordering and the corresponding supernodal nonzero
// structure of the implied factorization using the (Approximate) Minimum
// Degree reordering algorithm. Please see [ADD-96] for details.
//
// The input graph must be explicitly symmetric.
//
MinimumDegreeAnalysis MinimumDegree(
  const CoordinateGraph& graph, const MinimumDegreeControl& control);

} // namespace quotient

#include "quotient/minimum_degree-impl.hpp"

#endif // ifndef QUOTIENT_MINIMUM_DEGREE_H_
