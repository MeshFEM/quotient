/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_MACROS_H_
#define QUOTIENT_MACROS_H_

#ifdef _OPENMP
#include "omp.h"
#endif

#include <iostream>

// A guard for OpenMP pragmas so that builds which have not enabled OpenMP
// do not lead to compiler warnings.
#ifdef _OPENMP
#define OMP_PRAGMA(x) _Pragma(#x)
#else
#define OMP_PRAGMA(x)
#endif

// An attribute for routines which are known to not throw exceptions in
// release mode.
#ifdef QUOTIENT_DEBUG
#define QUOTIENT_NOEXCEPT
#else
#define QUOTIENT_NOEXCEPT noexcept
#endif

#ifdef QUOTIENT_DEBUG
#define QUOTIENT_ASSERT(assertion, message) \
  {                                         \
    if (!(assertion)) {                     \
      std::cerr << (message) << std::endl;  \
    }                                       \
  }
#else
#define QUOTIENT_ASSERT(condition, message)
#endif

// In most scenarios, it seems that the cost of the stronger hash is not
// worth the increased cost of maintaining the hash.
#ifdef QUOTIENT_STRONG_HASHES
#define QUOTIENT_HASH_COMBINE(hash, update) \
  hash ^= (update) + 0x9e3779b9 + (hash << 6) + (hash >> 2)
#else
#define QUOTIENT_HASH_COMBINE(hash, update) hash += update
#endif

#ifdef __GNUG__
#define QUOTIENT_UNUSED __attribute__((unused))
#elif defined(__clang__)
#define QUOTIENT_UNUSED __attribute__((unused))
#else
#define QUOTIENT_UNUSED
#endif  // ifdef __GNUG__

// A macro for faciliating the storage of two separate non-negative integers
// within a single signed integer. The mapping is an involution between
// [0, INT_MAX - 2] and -[2, INT_MAX].
//
// This choice allows -1 to be used as a special code.
#define SYMMETRIC_INDEX(index) (-((index) + 2))

#endif  // ifndef QUOTIENT_MACROS_H_
