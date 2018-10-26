/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_DEGREE_LISTS_IMPL_H_
#define QUOTIENT_DEGREE_LISTS_IMPL_H_

#include "quotient/config.hpp"
#include "quotient/degree_lists.hpp"

namespace quotient {

inline Int DegreeLists::FindMinimalIndex(bool demand_smallest_index) {
  while (degree_lower_bound < Int(degree_heads.size()) &&
      degree_heads[degree_lower_bound] == -1) {
    ++degree_lower_bound;
  }
  if (degree_lower_bound == Int(degree_heads.size())) {
    return -1;
  }

  Int index = degree_heads[degree_lower_bound];
  if (demand_smallest_index) {
    Int minimal_index = index;
    while (next_degree_member[index] != -1) {
      index = next_degree_member[index];
      minimal_index = std::min(minimal_index, index);
    }
    index = minimal_index;
  }
  return index;
}

inline void DegreeLists::RemoveDegree(Int index) {
  const Int degree = degrees[index];
  const Int last = last_degree_member[index];
  const Int next = next_degree_member[index];
  if (last == -1) {
    degree_heads[degree] = next;
  } else {
    // 'index' was not the head, so simply patch the connections.
    next_degree_member[last] = next;
  }
  if (next != -1) {
    last_degree_member[next] = last;
  }
}

inline void DegreeLists::AddDegree(Int index, Int degree) {
  const Int head = degree_heads[degree];
  degree_heads[degree] = index;
  last_degree_member[index] = -1;
  next_degree_member[index] = head;
  if (head != -1) {
    last_degree_member[head] = index;
  }
  degrees[index] = degree;

  // Update the minimal degree.
  degree_lower_bound = std::min(degree_lower_bound, degree);
}

inline void DegreeLists::UpdateDegree(Int index, Int degree) {
  if (degrees[index] == degree) {
    return;
  }
  RemoveDegree(index);
  AddDegree(index, degree);
}

} // namespace quotient

#endif // ifndef QUOTIENT_DEGREE_LISTS_IMPL_H_
