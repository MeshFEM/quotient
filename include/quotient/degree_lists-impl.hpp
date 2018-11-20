/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_DEGREE_LISTS_IMPL_H_
#define QUOTIENT_DEGREE_LISTS_IMPL_H_

#include "quotient/integers.hpp"
#include "quotient/macros.hpp"

#include "quotient/degree_lists.hpp"

namespace quotient {

inline Int DegreeLists::FindMinimalIndex(bool demand_smallest_index)
    QUOTIENT_NOEXCEPT {
  while (degree_lower_bound < Int(heads.size()) &&
         heads[degree_lower_bound] == -1) {
    ++degree_lower_bound;
  }
  QUOTIENT_ASSERT(degree_lower_bound != Int(heads.size()),
                  "Could not find a minimal degree.");

  Int index = heads[degree_lower_bound];
  if (demand_smallest_index) {
    Int minimal_index = index;
    while (next_member[index] != -1) {
      index = next_member[index];
      minimal_index = std::min(minimal_index, index);
    }
    index = minimal_index;
  }

#ifdef QUOTIENT_DEBUG
  Int degree = 0;
  while (heads[degree] == -1) {
    ++degree;
  }
  QUOTIENT_ASSERT(
      degree == degree_lower_bound,
      "True minimal degree was different from result from FindMinimalIndex.");
#endif

  RemoveDegree(index);

  return index;
}

inline void DegreeLists::RemoveDegree(Int index) QUOTIENT_NOEXCEPT {
  const Int degree = degrees[index];
  const Int last = last_member[index];
  const Int next = next_member[index];
  if (last == -1) {
    heads[degree] = next;
  } else {
    // 'index' was not the head, so simply patch the connections.
    next_member[last] = next;
  }
  if (next != -1) {
    last_member[next] = last;
  }
}

inline void DegreeLists::AddDegree(Int index, Int degree) QUOTIENT_NOEXCEPT {
  const Int head = heads[degree];
  heads[degree] = index;
  last_member[index] = -1;
  next_member[index] = head;
  if (head != -1) {
    last_member[head] = index;
  }
  degrees[index] = degree;

  // Update the minimal degree.
  degree_lower_bound = std::min(degree_lower_bound, degree);
}

inline void DegreeLists::UpdateDegree(Int index, Int degree) QUOTIENT_NOEXCEPT {
  if (degrees[index] == degree) {
    return;
  }
  RemoveDegree(index);
  AddDegree(index, degree);
}

}  // namespace quotient

#endif  // ifndef QUOTIENT_DEGREE_LISTS_IMPL_H_
