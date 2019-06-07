/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_IO_UTILS_IMPL_H_
#define QUOTIENT_IO_UTILS_IMPL_H_

#include "quotient/io_utils.hpp"

namespace quotient {

inline void ForestToDot(const std::string& filename,
                        const Buffer<Int>& parents) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Could not open " << filename << std::endl;
    return;
  }

  file << "digraph g{\n";
  for (std::size_t i = 0; i < parents.Size(); ++i) {
    if (parents[i] < 0) {
      continue;
    }
    std::ostringstream os;
    os << "  " << parents[i] << " -> " << i << ";\n";
    file << os.str();
  }
  file << "}\n";
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
  for (UInt i = 0; i < vec.size(); ++i) {
    os << vec[i] << " ";
  }
  return os;
}

template <typename T>
void Print(const std::vector<T>& vec, const std::string& msg,
           std::ostream& os) {
  os << msg << ": ";
  os << vec << std::endl;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Buffer<T>& vec) {
  for (UInt i = 0; i < vec.Size(); ++i) {
    os << vec[i] << " ";
  }
  return os;
}

template <typename T>
void Print(const Buffer<T>& vec, const std::string& msg, std::ostream& os) {
  os << msg << ": ";
  os << vec << std::endl;
}

}  // namespace quotient

#endif  // ifndef QUOTIENT_IO_UTILS_IMPL_H_
