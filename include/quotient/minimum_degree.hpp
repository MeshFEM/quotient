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

  // Whether nontrivial supernodes are allowed. It is highly recommended that
  // this remain true.
  bool allow_supernodes = true;

  // Whether aggressive element absorptions are allowed.
  bool aggressive_absorption = true;

  // Whether the list of pairs of aggressive element absorptions should be
  // returned in the MinimumDegreeAnalysis result of MinimumDegree.
  bool store_aggressive_absorptions = false;

  // Whether the list of pairs of variable merges should be returned in the
  // MinimumDegreeAnalysis result of MinimumDegree.
  bool store_variable_merges = false;

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
