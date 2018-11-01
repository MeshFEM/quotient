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
#include "quotient/minimum_degree_control.hpp"

namespace quotient {

// The result of running the MinimumDegree reordering algorithm. It contains
// the ordered list of eliminated principal vertices, the list of supernodes,
// and the supernodal nonzero structure of each principal column.
struct MinimumDegreeResult {
  // The recommended elimination order of the supernodes (with each supernode
  // represented by its principal member).
  std::vector<Int> elimination_order;

  // The list of supernodes. Entry 'i' is nonempty if and only if 'i' is the
  // principal member of the supernode, in which case 'supernodes[i]' contains
  // the sorted list of members of the supernode.
  std::vector<std::vector<Int>> supernodes;

  // The structures of the eliminated supernodes. Entry 'index' corresponds to
  // the structure of supernode 'elimination_order[index]'. This will only be
  // formed if 'MinimumDegreeControl.store_structures' was true.
  std::vector<std::vector<Int>> eliminated_structures;

  // The postordering of the assembly forest.
  std::vector<Int> postorder;

  // The number of aggressive absorptions that occurred.
  Int num_aggressive_absorptions;

  // An optional list (based on the value of
  // 'MinimumDegreeControl.store_pivot_element_list_sizes') of the lengths of
  // the element lists of each pivot.
  std::vector<Int> pivot_element_list_sizes;

  // The number of degree updates performed during the minimum-degree analysis.
  Int num_degree_updates = 0;

  // An optional count (based on the value of
  // 'MinimumDegreeControl.store_num_degree_updates_with_multiple_elements')
  // of the number of external degree updates that involved a variable with
  // more than two members in its element list.
  Int num_degree_updates_with_multiple_elements = -1;

  // The number of times that supervariables were falsely placed into the
  // same bucket.
  Int num_hash_bucket_collisions;

  // The number of times that supervariables falsely had the same hash value.
  Int num_hash_collisions;

  // The number of structural nonzeros in the Cholesky factor.
  Int num_cholesky_nonzeros = 0;

  // The number of floating-point operations in a standard Cholesky
  // factorization using this ordering.
  double num_cholesky_flops = 0;

  // An optional (based upon the value of 'MinimumDegreeControl.time_stages')
  // map from the stage names to the corresponding elapsed seconds.
  std::unordered_map<std::string, double> elapsed_seconds;

  // A trivial constructor.
  MinimumDegreeResult();

  // Returns the principal member of the largest supernode.
  Int LargestSupernode() const;

  // Returns the number of members of the largest supernode.
  Int LargestSupernodeSize() const;

  // Returns the number of strictly-lower nonzeros in the associated Cholesky
  // factor.
  Int NumStrictlyLowerCholeskyNonzeros() const;

  // Returns the fraction of pivots whose element list had more than two
  // members.
  double FractionOfPivotsWithMultipleElements() const;

  // Returns the fraction of degree updates whose corresponding variable had an
  // element list had more than two members.
  double FractionOfDegreeUpdatesWithMultipleElements() const;

  // Returns the permutation implied by the postordering of the assembly tree.
  std::vector<Int> Permutation() const;
};

// Returns a supernodal reordering and the corresponding supernodal nonzero
// structure of the implied factorization using the (Approximate) Minimum
// Degree reordering algorithm. Please see [ADD-96] for details.
//
// The input graph must be explicitly symmetric.
//
MinimumDegreeResult MinimumDegree(
  const CoordinateGraph& graph, const MinimumDegreeControl& control);

} // namespace quotient

#include "quotient/minimum_degree-impl.hpp"

#endif // ifndef QUOTIENT_MINIMUM_DEGREE_H_
