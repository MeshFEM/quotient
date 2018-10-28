/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_MACROS_H_
#define QUOTIENT_MACROS_H_

#include <cstdint>

namespace quotient {

#ifdef QUOTIENT_DEBUG
# define QUOTIENT_ASSERT(assertion, message) { \
  if (!(assertion)) { \
    std::cerr << (message) << std::endl; \
  } }
#else
# define QUOTIENT_ASSERT(condition, message)
#endif

// The datatype for signed indices.
typedef int64_t Int;

// The datatype used for unsigned indices.
typedef uint64_t UInt;

} // namespace quotient

#endif // ifndef QUOTIENT_MACROS_H_
