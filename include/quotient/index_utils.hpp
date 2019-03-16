/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_INDEX_UTILS_H_
#define QUOTIENT_INDEX_UTILS_H_

#include <vector>

#include "quotient/buffer.hpp"
#include "quotient/integers.hpp"

namespace quotient {

// Forms the inverse of a given permutation.
void InvertPermutation(const Buffer<Int>& permutation,
                       Buffer<Int>* inverse_permutation);

// Fills 'offsets' with a length 'num_indices + 1' array whose i'th index is
// the sum of the sizes whose indices are less than i.
void OffsetScan(const Buffer<Int>& sizes, Buffer<Int>* offsets);

// Builds the packed list of children of the nodes in a forest from the parent
// links. The children of node 'i' will be stored between indices
// 'child_offsets[i]' and 'child_offsets[i + 1]' of 'children'.
void ChildrenFromParents(const Buffer<Int>& parents, Buffer<Int>* children,
                         Buffer<Int>* child_offsets);

// Builds the packed list of children of the nodes in a given subsequence
// based upon their parents, which only exist at index 'i' if
// SYMMETRIC_INDEX(i) >= 0..
//
// As for 'ChildrenFromParents', the children of node 'i' will be stored
// between indices 'child_offsets[i]' and 'child_offsets[i + 1]' of 'children'.
void ChildrenFromParentSubsequence(const Buffer<Int>& symmetric_parents,
                                   const std::vector<Int>& node_subsequence,
                                   Buffer<Int>* children,
                                   Buffer<Int>* child_offsets);

}  // namespace quotient

#include "quotient/index_utils-impl.hpp"

#endif  // ifndef QUOTIENT_INDEX_UTILS_H_
