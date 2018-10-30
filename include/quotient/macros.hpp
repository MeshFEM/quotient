/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_MACROS_H_
#define QUOTIENT_MACROS_H_

namespace quotient {

#ifdef QUOTIENT_DEBUG
#  define QUOTIENT_ASSERT(assertion, message) { \
       if (!(assertion)) { \
         std::cerr << (message) << std::endl; \
       } }
#else
#  define QUOTIENT_ASSERT(condition, message)
#endif

// In most scenarios, it seems that the cost of the stronger hash is not
// worth the increased cost of maintaining the hash.
#ifdef QUOTIENT_STRONG_HASHES
#  define QUOTIENT_HASH_COMBINE(hash, update) \
       hash ^= (update) + 0x9e3779b9 + (hash << 6) + (hash >> 2)
#else
#  define QUOTIENT_HASH_COMBINE(hash, update) \
       hash += update
#endif

} // namespace quotient

#endif // ifndef QUOTIENT_MACROS_H_
