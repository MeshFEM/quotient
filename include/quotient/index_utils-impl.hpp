/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_INDEX_UTILS_IMPL_H_
#define QUOTIENT_INDEX_UTILS_IMPL_H_

#include "quotient/index_utils.hpp"

namespace quotient {

inline void InvertPermutation(const std::vector<Int>& permutation,
                              std::vector<Int>* inverse_permutation) {
  const Int num_indices = permutation.size();
  inverse_permutation->clear();
  inverse_permutation->resize(num_indices);
  for (Int index = 0; index < num_indices; ++index) {
    (*inverse_permutation)[permutation[index]] = index;
  }
}

inline void OffsetScan(const std::vector<Int>& sizes,
                       std::vector<Int>* offsets) {
  const Int num_indices = sizes.size();
  offsets->resize(num_indices + 1);

  Int offset = 0;
  for (Int index = 0; index < num_indices; ++index) {
    (*offsets)[index] = offset;
    offset += sizes[index];
  }
  (*offsets)[num_indices] = offset;
}

inline void ChildrenFromParents(const std::vector<Int>& parents,
                                std::vector<Int>* children,
                                std::vector<Int>* child_offsets) {
  const Int num_indices = parents.size();

  {
    // Compute the number of children of each node in the forest.
    std::vector<Int> num_children(num_indices, 0);
    for (Int index = 0; index < num_indices; ++index) {
      const Int parent = parents[index];
      if (parent >= 0) {
        ++num_children[parent];
      }
    }

    // Convert the number of children into the offsets to pack the children
    // into.
    OffsetScan(num_children, child_offsets);
  }

  // Pack the children into the offsets.
  children->resize(num_indices);
  auto offsets_copy = *child_offsets;
  for (Int index = 0; index < num_indices; ++index) {
    const Int parent = parents[index];
    if (parent >= 0) {
      (*children)[offsets_copy[parent]++] = index;
    }
  }
}

inline void ChildrenFromParentSubsequence(
    const std::vector<Int>& parents, const std::vector<Int>& node_subsequence,
    std::vector<Int>* children, std::vector<Int>* child_offsets) {
  const Int num_indices = parents.size();

  {
    // Compute the number of children of each node in the forest.
    std::vector<Int> num_children(num_indices, 0);
    for (const Int& index : node_subsequence) {
      const Int parent = parents[index];
      if (parent >= 0) {
        ++num_children[parent];
      }
    }

    // Convert the number of children into the offsets to pack the children
    // into.
    OffsetScan(num_children, child_offsets);
  }

  // Pack the children into the offsets.
  children->resize(num_indices);
  auto offsets_copy = *child_offsets;
  for (const Int& index : node_subsequence) {
    const Int parent = parents[index];
    if (parent >= 0) {
      (*children)[offsets_copy[parent]++] = index;
    }
  }
}

}  // namespace quotient

#endif  // ifndef QUOTIENT_INDEX_UTILS_IMPL_H_
