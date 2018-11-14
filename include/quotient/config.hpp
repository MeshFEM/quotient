/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_CONFIG_H_
#define QUOTIENT_CONFIG_H_

#ifdef _OPENMP
#include "omp.h"
#endif

#include <cstdint>

#ifdef _OPENMP
# define OMP_PRAGMA(x) _Pragma(#x)
#else
# define OMP_PRAGMA(x)
#endif

#ifdef QUOTIENT_DEBUG
# define QUOTIENT_NOEXCEPT
#else
// For now, we disable 'noexcept', even in release mode.
# define QUOTIENT_NOEXCEPT
#endif

namespace quotient {

// The datatype for signed indices.
#ifdef QUOTIENT_USE_64BIT
typedef int64_t Int;
#else
typedef int Int;
#endif

// The datatype used for unsigned indices.
#ifdef QUOTIENT_USE_64BIT
typedef uint64_t UInt;
#else
typedef unsigned UInt;
#endif

} // namespace quotient

#endif // ifndef QUOTIENT_CONFIG_H_
