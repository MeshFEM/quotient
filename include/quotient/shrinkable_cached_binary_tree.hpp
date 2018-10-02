/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_SHRINKABLE_CACHED_BINARY_TREE_H_
#define QUOTIENT_SHRINKABLE_CACHED_BINARY_TREE_H_

#include <iostream>
#include <vector>

#include "quotient/config.hpp"

namespace quotient {

inline Int PowerOfTwo(UInt exponent) {
  return UInt(1) << exponent;
}

inline UInt CeilLog2(UInt n) {
  UInt ceil_log2 = 0;
  while (n) {
    ++ceil_log2;
    n >>= UInt(1);
  }
  return ceil_log2;
}

// A data structure that, given an array of 'n' items of type 'T', supports:
//   * Random access lookups in O(1) time.
//   * Retrieval of an index with minimal value in O(1) time.
//   * Value modification in O(lg(n)) time.
//
// The associated storage mechanisms are:
//   * The list of 'n' values of type 'T'.
//   * The integer list of cached pointers to a minimal (reordered) child
//     leaf index.
//   * The permutation performed to make the kept indices contiguous.
//
// Nodes can be 'deleted' in O(lg(n)) time by swapping the 'deleted' index with
// the last active index, updating the two ancestor caches, updating the
// permutation array, and shrinking the variable denoting the number of kept
// elements.
//
// TODO(Jack Poulson): Provide support for overloading operator <.
//
template<typename T>
class ShrinkableCachedBinaryTree {
  // A trivial constructor.
  ShrinkableCachedBinaryTree();

  // A trivial destructor.
  ~ShrinkableCachedBinaryTree();

  // Initializes the cached tree based upon a given list of values.
  //
  // If there are 'n' members of 'values', this should require O(n) time.
  void Initialize(const std::vector<T>& values);

  // Retrieves the index and value of a minimal active entry.
  //
  // If there are 'num_active' active indices, then this should require
  // O(lg(num_active)) time to update the minimal value cache.
  std::pair<Int, T> MinimalEntry() const;

  // Retrieves the value of an index in O(1) time.
  const T& GetValue(Int index) const;

  // Modifies the value of an active index.
  //
  // If there are 'num_active' active indices, then this should require
  // O(lg(num_active)) time to update the minimal value cache.
  void SetValue(Int index, const T& value);

  // Updates the value of an active index.
  //
  // If there are 'num_active' active indices, then this should require
  // O(lg(num_active)) time to update the minimal value cache.
  void UpdateValue(Int index, const T& value);

  // Implicitly deletes an index.
  //
  // If there are 'num_active' active indices, then this should require
  // O(lg(num_active)) time to update the minimal value cache.
  void DisableIndex(Int index);

 private:
   // Returns the offset into 'comparison_tree_' of the given level.
   Int LevelOffset(Int level) const;

   // Returns the number of entries in 'comparison_tree_' for the given level.
   Int LevelSize(Int level) const;

   // Returns the flattened tree index of a minimal entry in the child's
   // subtree.
   Int ChildMinimalTreeIndex(Int level, Int index, bool right_child) const;

   // Updates the entry of 'comparison_tree_' corresponding to the given level
   // and index within said level using its two children.
   void UpdateComparisonUsingChildren(Int level, Int index);

   // Returns a pair of the tree level and index of the parent node of the
   // given leaf node.
   std::pair<Int, Int> ParentLevelAndIndex(Int tree_index) const;

   // An ordered list of values of type 'T'.
   std::vector<T> values_;

   // The number of members of 'values_' that have not been implicitly deleted.
   Int num_active_;

   // A permutation over [0, values_.size()) such that the 'active' values
   // map to [0, num_active_).
   std::vector<Int> perm_; 

   // The inverse of 'perm_': a map over [0, values_.size()) that maps
   // [0, num_active_) into the original positions in [0, values_.size()).
   std::vector<Int> inverse_perm_;

   // An array of length 'num_active_ - 1' that stores the packed indices of
   // the minimal (permuted) active index of the minimal descendant.
   //
   // For example, if there are four active indices, say, with values:
   //   [7, 2, 5, 3],
   // then 'comparison_tree_' will store a packed version of the binary tree:
   //
   //           1
   //          / \
   //         -   -
   //       1       3
   //      / \     / \
   //     0   1   2   3
   //
   // The packing would be the serialization of [[1], [1, 3]], where the root
   // node (level 0) is stored first, then the nodes at level 1, etc. And, in
   // cases where there are not a power of two number of nodes, the bottom-most
   // level is truncated so that there are exactly 'num_active_ - 1' nodes in
   // the tree. For example, if there are five active indices, say, with values:
   //   [4, 8, 0, 5, 2],
   // then we would have the binary tree:
   //
   //               2
   //              / \
   //            --   --
   //          2         4         
   //         / \       / \
   //        0   2     3   4
   //       / \
   //      0   1
   //
   // The packing would be the serialization of [[2], [2, 4], [0]].
   //
   // Programatically, the j'th level of the tree is stored in indices:
   //   [2^j - 1, ..., max(2^{j + 1} - 2, num_active_ - 1)].
   std::vector<Int> comparison_tree_;

   // The number of levels of comparison metadata.
   Int num_comparison_levels_;
};

template<typename T>
ShrinkableBinaryTree<T>::ShrinkableBinaryTree() { }

template<typename T>
ShrinkableBinaryTree<T>::~ShrinkableBinaryTree() { }

template<typename T>
Int ShrinkableBinaryTree<T>::LevelOffset(Int level) const {
  return PowerOfTwo(level) - 1;
}

template<typename T>
Int ShrinkableBinaryTree<T>::LevelSize(Int level) const {
#ifdef QUOTIENT_DEBUG
  if (level < 0 || level >= num_comparison_levels_) {
    std::cerr << "Requested level size of invalid level." << std::endl;
    return 0;
  }
#endif
  if (level == num_comparison_levels_ - 1) {
    return num_active_ - PowerOfTwo(num_comparison_levels_ - 1);
  }
  return PowerOfTwo(level);
}

template<typename T>
void ShrinkableBinaryTree<T>::ChildMinimalTreeIndex(
    Int level, Int index, bool right_child) const {
#ifdef QUOTIENT_DEBUG
  if (level < 0 || level >= num_comparison_levels_) {
    std::cerr << "Requested level size of invalid level." << std::endl;
    return;
  }
  if (index < 0 || index >= LevelSize(level)) {
    std::cerr << "Requested invalid index of level." << std::endl;
    return;
  }
#endif
  const Int child_index = right_child ? 2 * index + 1 : 2 * index;
  if (level == num_comparison_levels_ - 1) {
    return child_index;
  }
  if (child_index >= LevelSize(level + 1)) {
    // This should only be possible when level is num_comparison_levels_ - 2.
    return 2 * last_level_size + (child_index - last_level_size);
  }
  return comparison_tree_[LevelOffset(level + 1) + child_index];
}

template<typename T>
bool ShrinkableBinaryTree<T>::UpdateComparisonUsingChildren(
    Int level, Int index) {
  const Int left_tree_index = ChildMinimalTreeIndex(
      level, index, false /* right_child */);
  const Int right_tree_index = ChildMinimalTreeIndex(
      level, index, true /* right_child */);
  const Int level_offset = LevelOffset(level);
  const bool left_is_small =
      values_[inverse_perm_[left_tree_index]] <=
      values_[inverse_perm_[right_tree_index]];
  const Int new_index = left_is_small ? left_tree_index : right_tree_index;
  const Int old_index = comparison_tree_[level_offset + index];
  comparison_tree_[level_offset + index] = new_index;
  return new_index == old_index;
}

template<typename T>
void ShrinkableCachedBinaryTree<T>::Initialize(const std::vector<T>& values) {
  values_ = values;
  num_active_ = values.size();
  num_comparison_levels_ = CeilLog2(num_active_);
  if (num_active_ == 0) {
    return;
  }

  // Create an identity permutation
  perm_.resize(values.size());
  std::iota(perm_.begin(), perm_.end(), 0);
  inverse_perm_ = perm_;

  // Compute the comparison metadata from the leaves up (in linear time).
  comparison_tree_.resize(num_active_ - 1);
  for (Int level = num_comparison_levels - 1; level >= 0; --level) {
    const Int level_size = LevelSize(level);
    for (Int index = 0; index < level_size; ++index) {
      UpdateComparisonUsingChildren(level, index);
    }
  }
}

template<typename T>
std::pair<Int, T> ShrinkableCachedBinaryTree<T>::MinimalEntry() const {
  std::pair<Int, T> entry;
  entry.first = inverse_perm_[comparison_tree_[0]];
  entry.second = values_[entry.first];
  return entry;
}

template<typename T>
std::pair<Int, Int> ShrinkableCachedBinaryTree<T>::ParentLevelAndIndex(
    Int tree_index) const {
  std::pair<Int, Int> tree_pos;
  const Int last_level_size = LevelSize(num_comparison_levels_ - 1);
  if (tree_index < 2 * last_level_size) {
    tree_pos.first = num_comparison_levels_ - 1;
    tree_pos.second = tree_index / 2;
  } else {
    tree_pos.first = num_comparison_levels_ - 2;
    tree_pos.second = tree_index - 2 * last_level_size;
  }
  return tree_pos;
}

template<typename T>
void ShrinkableCachedBinaryTree<T>::PropagateComparisons(Int index) {
  const Int tree_index = perm_[index];
  std::pair<Int, Int> tree_pos = ParentLevelAndIndex(tree_index);
  bool done = UpdateComparisonUsingChildren(tree_pos.first, tree_pos.second);
  while (!done && tree_pos.first != 0) {
    --tree_pos.first;
    tree_pos.second /= 2;
    done = UpdateComparisonUsingChildren(tree_pos.first, tree_pos.second);
  }
}

template<typename T>
const T& ShrinkableCachedBinaryTree<T>::GetValue(Int index) const {
  return values_[index];
}

template<typename T>
void ShrinkableCachedBinaryTree<T>::SetValue(Int index, const T& value) {
  if (values_[index] == value) {
    return;
  }
  values_[index] = value;
  PropagateComparisons(index);
}

template<typename T>
void ShrinkableCachedBinaryTree<T>::UpdateValue(Int index, const T& value) {
  if (value == T(0)) {
    return;
  }
  values_[index] += value;
  PropagateComparisons(index);
}

template<typename T>
void ShrinkableCachedBinaryTree<T>::DisableIndex(Int index) {
  const Int tree_index = perm_[index];
#ifdef QUOTIENT_DEBUG
  if (tree_index >= num_active_) {
    std::cerr << "Index was already disabled." << std::endl;
    return;
  }
#endif
  const Int index0 = inverse_perm_[tree_index];
  const Int index1 = inverse_perm_[num_active_ - 1];
  if (tree_index != num_active_ - 1) {
    std::swap(inverse_perm_[tree_index], inverse_perm_[num_active_ - 1]);
    std::swap(perm_[index0], perm_[index1]);
  }

  --num_active_;
  num_comparison_levels_ = CeilLog2(num_active_);
  comparison_tree_.pop_back();

  PropagateComparisons(tree_index);
  PropagateComparisons(num_active - 1);
}

} // namespace quotient

#endif // ifndef QUOTIENT_SHRINKABLE_CACHED_BINARY_TREE_H_
