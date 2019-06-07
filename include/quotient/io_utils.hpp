/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_IO_UTILS_H_
#define QUOTIENT_IO_UTILS_H_

#include <ostream>
#include <string>
#include <vector>

#include "quotient/buffer.hpp"
#include "quotient/integers.hpp"

namespace quotient {

// Writes a dot file (usually ".gv") for the forest implied by the parents.
// One can subsequently generate a PNG of the forest using:
//   dot -Tpng filename -o output.png
// But beware that the call to dot might take 15 minutes or more.
void ForestToDot(const std::string& filename, const Buffer<Int>& parents);

// Pretty-prints an std::vector<T>.
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec);
template <typename T>
void Print(const std::vector<T>& vec, const std::string& msg, std::ostream& os);

// Pretty-prints a Buffer<T>.
template <typename T>
std::ostream& operator<<(std::ostream& os, const Buffer<T>& vec);
template <typename T>
void Print(const Buffer<T>& vec, const std::string& msg, std::ostream& os);

}  // namespace quotient

#include "quotient/io_utils-impl.hpp"

#endif  // ifndef QUOTIENT_IO_UTILS_H_
