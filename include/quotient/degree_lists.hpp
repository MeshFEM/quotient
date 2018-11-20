/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_DEGREE_LISTS_H_
#define QUOTIENT_DEGREE_LISTS_H_

#include <vector>

#include "quotient/integers.hpp"
#include "quotient/macros.hpp"

namespace quotient {

// An array of  doubly-linked list for tracking node degrees in a manner which
// makes it easy to find a minimal degree variable.
struct DegreeLists {
  // A lower-bound on the minimum degree of all degree list members.
  Int degree_lower_bound = 0;

  // A list of the external degrees associated with each index.
  std::vector<Int> degrees;

  // The indices of the heads of the lists for each degree. If no index exists
  // for 'degree', then 'heads[degree]' is -1.
  std::vector<Int> heads;

  // Index 'i' will provide the index of the member that occurs directly after
  // 'i' in the degree list containing index 'i' (if such a member does not
  // exist, the value will be -1).
  std::vector<Int> next_member;

  // Index 'i' will provide the index of the member that occurs directly before
  // 'i' in the degree list containing index 'i' (if such a member does not
  // exist, the value will be -1).
  std::vector<Int> last_member;

  // Successively increases degree_lower_bound until a degree list member is
  // found (and an index of minimal degree is returned).
  // If 'demand_smallest_index' is true, then the member of the linked list with
  // the smallest index is returned.
  Int FindMinimalIndex(bool demand_smallest_index) QUOTIENT_NOEXCEPT;

  // Removes the current degree of occurrence with the given index.
  void RemoveDegree(Int index) QUOTIENT_NOEXCEPT;

  // Adds in an occurrence of the specified index and degree.
  void AddDegree(Int index, Int degree) QUOTIENT_NOEXCEPT;

  // If the old degree was the same as the current degree, then this routine
  // removes the old occurrence and adds in the new occurrence. Otherwise,
  // it is a no-op.
  void UpdateDegree(Int index, Int degree) QUOTIENT_NOEXCEPT;
};

}  // namespace quotient

#include "quotient/degree_lists-impl.hpp"

#endif  // ifndef QUOTIENT_DEGREE_LISTS_H_
