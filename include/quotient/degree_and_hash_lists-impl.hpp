/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_DEGREE_AND_HASH_LISTS_IMPL_H_
#define QUOTIENT_DEGREE_AND_HASH_LISTS_IMPL_H_

#include "quotient/integers.hpp"
#include "quotient/macros.hpp"

#include "quotient/degree_and_hash_lists.hpp"

namespace quotient {

inline Int DegreeAndHashLists::FindMinimalIndex(bool demand_smallest_index)
    QUOTIENT_NOEXCEPT {
  while (degree_lower_bound < Int(heads.Size()) &&
         heads[degree_lower_bound] == -1) {
    ++degree_lower_bound;
  }
  QUOTIENT_ASSERT(degree_lower_bound != Int(heads.Size()),
                  "Could not find a minimal degree.");

#ifdef QUOTIENT_DEBUG
  Int degree = 0;
  while (heads[degree] == -1) {
    ++degree;
  }
  QUOTIENT_ASSERT(
      degree == degree_lower_bound,
      "True minimal degree was different from result from FindMinimalIndex.");
#endif

  Int index = heads[degree_lower_bound];
  if (demand_smallest_index) {
    Int minimal_index = index;
    while (next_member[index] != -1) {
      index = next_member[index];
      minimal_index = std::min(minimal_index, index);
    }
    index = minimal_index;
    RemoveDegree(index);
  } else {
    RemoveHeadDegree(index, degree_lower_bound);
  }

  return index;
}

inline void DegreeAndHashLists::RemoveHeadDegree(Int index,
                                                 Int degree) QUOTIENT_NOEXCEPT {
  const Int next = next_member[index];
  QUOTIENT_ASSERT(last_member[index] == -1, "Falsely assumed head index.");
  heads[degree] = next;
  if (next != -1) {
    last_member[next] = -1;
  }
}

inline void DegreeAndHashLists::RemoveDegree(Int index) QUOTIENT_NOEXCEPT {
  QUOTIENT_ASSERT(index >= 0 && index < Int(heads.Size()),
                  "Index of " + std::to_string(index) + " was out of bounds.");
  const Int degree = degrees[index];
  QUOTIENT_ASSERT(
      degree >= 0 && degree < Int(heads.Size()),
      "Degree of " + std::to_string(degree) + " was out of bounds.");
  const Int last = last_member[index];
  const Int next = next_member[index];
  if (last == -1) {
    heads[degree] = next;
  } else {
    QUOTIENT_ASSERT(
        last >= 0 && last < Int(heads.Size()),
        "Last index of " + std::to_string(last) + " was out of bounds.");
    // 'index' was not the head, so simply patch the connections.
    next_member[last] = next;
  }
  if (next != -1) {
    QUOTIENT_ASSERT(
        next >= 0 && next < Int(heads.Size()),
        "Next index of " + std::to_string(next) + " was out of bounds.");
    last_member[next] = last;
  }
}

inline void DegreeAndHashLists::AddDegree(Int index,
                                          Int degree) QUOTIENT_NOEXCEPT {
  QUOTIENT_ASSERT(index >= 0 && index < Int(heads.Size()),
                  "Index of " + std::to_string(index) + " was out of bounds.");
  QUOTIENT_ASSERT(
      degree >= 0 && degree < Int(heads.Size()),
      "Degree of " + std::to_string(degree) + " was out of bounds.");
  const Int head = heads[degree];
  QUOTIENT_ASSERT(head != index, "Index matched preexisting head.");
  if (head != -1) {
    QUOTIENT_ASSERT(degrees[head] == degree, "Invalid head degree for degree " +
                                                 std::to_string(degree) +
                                                 " of " + std::to_string(head));
    last_member[head] = index;
  }
  heads[degree] = index;
  last_member[index] = -1;
  next_member[index] = head;
  degrees[index] = degree;

  // Update the minimal degree.
  degree_lower_bound = std::min(degree_lower_bound, degree);
}

inline void DegreeAndHashLists::UpdateDegree(Int index,
                                             Int degree) QUOTIENT_NOEXCEPT {
  if (degrees[index] == degree) {
    return;
  }
  RemoveDegree(index);
  AddDegree(index, degree);
}

inline Int DegreeAndHashLists::Hash(Int index) const QUOTIENT_NOEXCEPT {
  // The hash lists are only single-linked (using next_member), so we use the
  // back-link buffer (last_member) over the set of indices reserved for
  // hash lists (the principal members of the pivot structure).
  QUOTIENT_ASSERT(index >= 0 && index < Int(heads.Size()),
                  "Index of " + std::to_string(index) + " was out of bounds.");
  return last_member[index];
}

inline void DegreeAndHashLists::SetHash(Int index,
                                        std::size_t hash) QUOTIENT_NOEXCEPT {
  // The hash lists are only single-linked (using next_member), so we can
  // use the back-link buffer (last_member) over the set of indices reserved
  // for hash lists (the principal members of the pivot structure).
  QUOTIENT_ASSERT(index >= 0 && index < Int(heads.Size()),
                  "Index of " + std::to_string(index) + " was out of bounds.");
  last_member[index] = hash;
}

inline void DegreeAndHashLists::AddHash(Int index, std::size_t hash,
                                        Int bucket) QUOTIENT_NOEXCEPT {
  QUOTIENT_ASSERT(index >= 0 && index < Int(heads.Size()),
                  "Index of " + std::to_string(index) + " was out of bounds.");
  QUOTIENT_ASSERT(bucket >= 0,
                  "Bucket was negative: " + std::to_string(bucket));
  QUOTIENT_ASSERT(bucket < Int(heads.Size()),
                  "Bucket was out-of-scope. bucket: " + std::to_string(bucket) +
                      ", heads.Size(): " + std::to_string(heads.Size()));
  SetHash(index, hash);
  if (heads[bucket] >= 0) {
    // There is a non-empty degree list of degree 'bucket'.
    const Int head = last_member[heads[bucket]];
    QUOTIENT_ASSERT(index != head, "Tried to set next_member[" +
                                       std::to_string(index) +
                                       "] := " + std::to_string(index));
    next_member[index] = head;
    last_member[heads[bucket]] = index;
  } else {
    const Int head = SYMMETRIC_INDEX(heads[bucket]);
    QUOTIENT_ASSERT(index != head, "Tried to set next_member[" +
                                       std::to_string(index) +
                                       "] := " + std::to_string(index));
    next_member[index] = head;
    heads[bucket] = SYMMETRIC_INDEX(index);
  }
}

inline void DegreeAndHashLists::ClearHashBucket(Int bucket) QUOTIENT_NOEXCEPT {
  // Because we only traverse the next_member chain from the head of a
  // bucket until reaching -1, if we simply set the head of the bucket to -1,
  // then the bucket is effectively cleared.
  if (heads[bucket] >= 0) {
    // There is a non-empty degree list of degree 'bucket'.
    last_member[heads[bucket]] = -1;
  } else {
    heads[bucket] = -1;
  }
}

inline Int DegreeAndHashLists::HashBucketHead(Int bucket) const
    QUOTIENT_NOEXCEPT {
  QUOTIENT_ASSERT(bucket >= 0,
                  "Bucket was negative: " + std::to_string(bucket));
  if (heads[bucket] >= 0) {
    // There is a non-empty degree list of degree 'bucket'.
    return last_member[heads[bucket]];
  } else {
    return SYMMETRIC_INDEX(heads[bucket]);
  }
}

inline void DegreeAndHashLists::SetHashBucketHead(Int bucket,
                                                  Int head) QUOTIENT_NOEXCEPT {
  QUOTIENT_ASSERT(bucket >= 0,
                  "Bucket was negative: " + std::to_string(bucket));
  if (heads[bucket] >= 0) {
    // There is a non-empty degree list of degree 'bucket'.
    last_member[heads[bucket]] = head;
  } else {
    heads[bucket] = SYMMETRIC_INDEX(head);
  }
}

}  // namespace quotient

#endif  // ifndef QUOTIENT_DEGREE_AND_HASH_LISTS_IMPL_H_
