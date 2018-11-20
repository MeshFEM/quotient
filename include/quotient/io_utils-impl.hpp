/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_IO_UTILS_IMPL_H_
#define QUOTIENT_IO_UTILS_IMPL_H_

#include <ostream>
#include <vector>

#include "quotient/integers.hpp"
#include "quotient/io_utils.hpp"

namespace quotient {

template<typename T>
void PrintVector(
    const std::vector<T>& vec, const std::string& msg, std::ostream& os) {
  os << msg << ": ";
  for (UInt i = 0; i < vec.size(); ++i) {
    os << vec[i] << " ";
  }
  os << "\n";
}

} // namespace quotient

#endif // ifndef QUOTIENT_IO_UTILS_IMPL_H_
