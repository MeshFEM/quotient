/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_HASH_LISTS_H_
#define QUOTIENT_HASH_LISTS_H_

#include <vector>

#include "quotient/integers.hpp"
#include "quotient/macros.hpp"

namespace quotient {

// An array of singly-linked list for tracking variables with colliding hash
// bucket values.
struct HashLists {
  // A list of the hash buckets associated with each index.
  std::vector<Int> buckets;

  // A list of the hashes associated with each index.
  std::vector<std::size_t> hashes;

  // The indices of the heads of the lists for each bucket. If no
  // index exists for 'bucket', then 'heads[bucket]' is -1.
  std::vector<Int> heads;

  // Index 'i' will provide the index of the member that occurs directly after
  // 'i' in the list containing index 'i' (if such a member does not exist,
  // the value will be -1).
  std::vector<Int> next_member;

  // Removes the given hash bucket.
  void ClearBucket(Int bucket) QUOTIENT_NOEXCEPT;

  // Adds in an occurrence of the specified index and hash.
  void AddHash(Int index, std::size_t hash, Int bucket) QUOTIENT_NOEXCEPT;
};

} // namespace quotient

#include "quotient/hash_lists-impl.hpp"

#endif // ifndef QUOTIENT_HASH_LISTS_H_
