/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_MINIMUM_DEGREE_CONTROL_H_
#define QUOTIENT_MINIMUM_DEGREE_CONTROL_H_

#include "quotient/config.hpp"
#include "quotient/external_degree_type.hpp"

namespace quotient {

// A choice of hash function for mapping a variable to an std::size_t.
enum VariableHashType {
  // Due to Ashcraft and described in [ADD-96].
  kAshcraftVariableHash,

  // A custom hash function meant to avoid collisions by removing the
  // modular arithmetic (allowing the full range of std::size_t) and
  // multiplying each index update by its position in the adjaency or
  // element list.
  kBasicVariableHash,
};

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
  bool aggressive_absorption = false;

  // Whether the entire degree list will be traversed in order to ensure that
  // the member of the minimal degree list with smallest index is chosen.
  bool force_minimal_pivot_indices = false;

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

} // namespace quotient

#endif // ifndef QUOTIENT_MINIMUM_DEGREE_CONTROL_H_
