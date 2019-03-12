/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_HASH_LISTS_IMPL_H_
#define QUOTIENT_HASH_LISTS_IMPL_H_

#include "quotient/integers.hpp"
#include "quotient/macros.hpp"

#include "quotient/hash_lists.hpp"

namespace quotient {

inline void HashLists::AddHash(Int index, std::size_t hash,
                               Int bucket) QUOTIENT_NOEXCEPT {
  hashes[index] = hash;

  const Int head = heads[bucket];
  next_member[index] = head;
  heads[bucket] = index;
}

inline void HashLists::ClearBucket(Int bucket) QUOTIENT_NOEXCEPT {
  // Because we only traverse the next_member chain from the head of a
  // bucket until reaching -1, if we simply set the head of the bucket to -1,
  // then the bucket is effectively cleared.
  heads[bucket] = -1;
}

}  // namespace quotient

#endif  // ifndef QUOTIENT_HASH_LISTS_IMPL_H_
