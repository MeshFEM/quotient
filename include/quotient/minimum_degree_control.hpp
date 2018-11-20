/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_MINIMUM_DEGREE_CONTROL_H_
#define QUOTIENT_MINIMUM_DEGREE_CONTROL_H_

#include "quotient/degree_type.hpp"
#include "quotient/integers.hpp"

namespace quotient {

// A data structure for controlling the MinimumDegree reordering routine.
struct MinimumDegreeControl {
  // The type of approximation to use for the degree estimates.
  DegreeType degree_type = kAmestoyDegree;

  // Whether aggressive element absorptions are allowed.
  bool aggressive_absorption = true;

  // If a row has at least
  //
  //   max(min_dense_threshold,
  //       dense_sqrt_multiple * sqrt(num_original_vertices)),
  //
  // nonzeros away from the diagonal, then it will be treated as dense and
  // placed at the end of the ordering.
  Int min_dense_threshold = 16;
  float dense_sqrt_multiple = 10.f;

  // Whether the entire degree list will be traversed in order to ensure that
  // the member of the minimal degree list with smallest index is chosen.
  bool force_minimal_pivot_indices = false;

  // Whether a list should be stored of the lengths of the element lists of
  // the pivots.
  bool store_pivot_element_list_sizes = false;

  // Whether or not to store the count of the number of degree updates which
  // involved more than two (and, separately, how many less than or equal to
  // two) elements in the variable's element list.
  bool store_num_degree_updates_with_multiple_elements = false;

  // Whether nontrivial supernodes are allowed. It is highly recommended that
  // this remain true.
  bool allow_supernodes = true;

  // Return the explicit list of supernodes (in the original ordering)?
  bool store_supernodes = false;

  // Whether or not to create and store the nonzeros structures of each pivot.
  bool store_structures = false;
};

}  // namespace quotient

#endif  // ifndef QUOTIENT_MINIMUM_DEGREE_CONTROL_H_
