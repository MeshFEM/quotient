/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_RANDOM_ACCESS_HEAP_H_
#define QUOTIENT_RANDOM_ACCESS_HEAP_H_

#include <iostream>
#include <vector>

#include "quotient/config.hpp"

namespace quotient {

// A data structure that, given an array of 'n' items of type 'T', supports:
//   * Random access lookups in O(1) time.
//   * Retrieval of the left-most index with minimal value in O(1) time.
//   * Value modification in O(lg(n)) time.
//
// The associated storage mechanisms are:
//   * The list of 'n' values of type 'T'.
//   * The integer list of cached pointers to a minimal (reordered) child
//     leaf index.
//   * A boolean list of whether or not indices are still valid.
//
// Indices can be 'disabled' in O(lg(n)) time by setting the valid flag to
// false and then propagating the comparison changes up the tree.
//
// The high-level inspiration for the data structure is:
//
//   https://dxr.mozilla.org/mozilla-central/rev/1b822c694681/dom/canvas/WebGLElementArrayCache.cpp,
//
// which was written by Benoit Jacob to aid in preventing malicious OpenGL
// indexing attacks.
//
// TODO(Jack Poulson): Add support for a PushBack(const T& value) method and
// a DeleteIndex(Int index) method.
//
// TODO(Jack Poulson): Provide support for overloading operator <.
//
// TODO(Jack Poulson): Add support for aggolmerating the first few levels of
// the tree for performance reasons.
//
template<typename T>
class RandomAccessHeap {
 public:
  // A trivial constructor.
  RandomAccessHeap();

  // A trivial destructor.
  ~RandomAccessHeap();

  // Initializes the cached tree based upon a given list of values.
  //
  // If there are 'n' members of 'values', this should require O(n) time.
  void Reset(const std::vector<T>& values);

  // Retrieves the index and value of a minimal active entry.
  //
  // If there are 'num_active' active indices, then this should require
  // O(lg(num_active)) time to update the minimal value cache.
  std::pair<Int, T> MinimalEntry() const;

  // Retrieves a reference to the underlying vector of booleans denoting if
  // each value is currently valid.
  const std::vector<bool>& ValidValues() const;

  // Retrieves the validity of a given index.
  bool ValidValue(Int index) const;

  // Retrieves a reference to the underlying vector of values.
  const std::vector<T>& Values() const;

  // Retrieves the value of an index in O(1) time.
  const T& Value(Int index) const;

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

  // Returns 2^exponent.
  static Int PowerOfTwo(UInt exponent);

  // Returns ceil(log2(n)).
  static UInt CeilLog2(UInt n);

 private:
   // Returns the offset into 'comparison_tree_' of the given level.
   Int LevelOffset(Int level) const;

   // Returns the number of entries in 'comparison_tree_' for the given level.
   Int LevelSize(Int level) const;

   // Returns the flattened tree index of a minimal entry in the child's
   // subtree.
   Int ChildMinimalIndex(Int level, Int level_index, bool right_child) const;

   // Updates all relevant ancestral metadata for the leaf node with given
   // index.
   void PropagateComparisons(Int index);

   // In order to guarantee that the left-most (in the original ordering)
   // minimal entry is chosen in the case of a tie, we need to incorporate the
   // permutation into the tie-breaking decision.
   //
   // But we also avoid allowing an invalid index from being chosen (unless
   // both indices are invalid).
   bool UseLeftIndex(Int left_index, Int right_index) const;

   // Updates the entry of 'comparison_tree_' corresponding to the given level
   // and index within said level using its two children.
   void UpdateComparisonUsingChildren(Int level, Int level_index);

   // Returns true if the comparison metadata at the given tree position is
   // valid.
   bool ComparisonIsValid(Int level, Int level_index) const;

   // Returns true if the entire tree is valid.
   bool TreeIsValid() const;

   // Returns a pair of the tree level and index of the parent node of the
   // given leaf node.
   std::pair<Int, Int> ParentLevelAndIndex(Int index) const;

   // An ordered list of values of type 'T'.
   std::vector<T> values_;

   // An ordered list of booleans indicating if each member of 'values_' is
   // an acceptable minimum.
   std::vector<bool> valid_values_;

   /* An array of length 'num_values - 1' that stores the packed indices of
      the minimal (permuted) active index of the minimal descendant.

      For example, if there are four valid indices, say, with values:
        [7, 2, 5, 3],
      then 'comparison_tree_' will store a packed version of the binary tree:
   
                1
               / \
              -   -
            1       3
           / \     / \
          0   1   2   3
     
      The packing would be the serialization of [[1], [1, 3]], where the root
      node (level 0) is stored first, then the nodes at level 1, etc. And, in
      cases where there are not a power of two number of nodes, the bottom-most
      level is truncated so that there are exactly 'num_values - 1' nodes in
      the tree. For example, if there are five valid indices, say, with values:
        [4, 8, 0, 5, 2],
      then we would have the binary tree:
   
                    2
                   / \
                 --   --
               2         4         
              / \       / \
             0   2     3   4
            / \
           0   1
     
      The packing would be the serialization of [[2], [2, 4], [0]].
     
      Programatically, the j'th level of the tree is stored in indices:
        [2^j - 1, ..., max(2^{j + 1} - 2, num_values - 1)].
   */
   std::vector<Int> comparison_tree_;

   // The number of levels of comparison metadata.
   Int num_comparison_levels_;
};

template<typename T>
RandomAccessHeap<T>::RandomAccessHeap()
: num_comparison_levels_(0) { }

template<typename T>
RandomAccessHeap<T>::~RandomAccessHeap() { }

template<typename T>
Int RandomAccessHeap<T>::LevelOffset(Int level) const {
  return PowerOfTwo(level) - 1;
}

template<typename T>
Int RandomAccessHeap<T>::LevelSize(Int level) const {
#ifdef QUOTIENT_DEBUG
  if (level < 0 || level >= num_comparison_levels_) {
    std::cerr << "Requested level size of invalid level." << std::endl;
    return 0;
  }
#endif
  if (level == num_comparison_levels_ - 1) {
    return values_.size() - PowerOfTwo(num_comparison_levels_ - 1);
  }
  return PowerOfTwo(level);
}

template<typename T>
Int RandomAccessHeap<T>::ChildMinimalIndex(
    Int level, Int level_index, bool right_child) const {
#ifdef QUOTIENT_DEBUG
  if (level < 0 || level >= num_comparison_levels_) {
    std::cerr << "Requested level size of invalid level." << std::endl;
    return -1;
  }
  if (level_index < 0 || level_index >= LevelSize(level)) {
    std::cerr << "Requested invalid index of level." << std::endl;
    return -1;
  }
#endif
  const Int child_index = right_child ? 2 * level_index + 1 : 2 * level_index;
  if (level == num_comparison_levels_ - 1) {
    return child_index;
  }
  if (child_index >= LevelSize(level + 1)) {
    // This should only be possible when level is num_comparison_levels_ - 2.
#ifdef QUOTIENT_DEBUG
    if (level != num_comparison_levels_ - 2) {
      std::cerr << "Impossible level comparison." << std::endl;
    }
#endif
    const Int last_level_size = LevelSize(level + 1);
    return 2 * last_level_size + (child_index - last_level_size);
  }
  return comparison_tree_[LevelOffset(level + 1) + child_index];
}

template<typename T>
bool RandomAccessHeap<T>::UseLeftIndex(Int left_index, Int right_index) const {
  bool use_left_index;
  if (valid_values_[left_index] && valid_values_[right_index]) {
    use_left_index = values_[left_index] <= values_[right_index];
  } else if (valid_values_[left_index]) {
    use_left_index = true;
  } else if (valid_values_[right_index]) {
    use_left_index = false;
  } else {
    // Fall back on the left index.
    use_left_index = true;
  }
  return use_left_index;
}

template<typename T>
void RandomAccessHeap<T>::UpdateComparisonUsingChildren(
    Int level, Int level_index) {
  const Int left_index = ChildMinimalIndex(
      level, level_index, false /* right_child */);
  const Int right_index = ChildMinimalIndex(
      level, level_index, true /* right_child */);

  const bool use_left_index = UseLeftIndex(left_index, right_index);
  const Int new_index = use_left_index ? left_index : right_index;

  const Int level_offset = LevelOffset(level);
  comparison_tree_[level_offset + level_index] = new_index;
}

template<typename T>
bool RandomAccessHeap<T>::ComparisonIsValid(Int level, Int level_index)  const {
  const Int left_index = ChildMinimalIndex(
      level, level_index, false /* right_child */);
  const Int right_index = ChildMinimalIndex(
      level, level_index, true /* right_child */);

  const bool use_left_index = UseLeftIndex(left_index, right_index);
  const Int valid_index = use_left_index ? left_index : right_index;

  const Int level_offset = LevelOffset(level);
  return comparison_tree_[level_offset + level_index] == valid_index;
}

template<typename T>
bool RandomAccessHeap<T>::TreeIsValid() const {
  for (Int level = num_comparison_levels_ - 1; level >= 0; --level) {
    const Int level_size = LevelSize(level);
    for (Int level_index = 0; level_index < level_size; ++level_index) {
      if (!ComparisonIsValid(level, level_index)) {
#ifdef QUOTIENT_DEBUG
        std::cerr << "Comparison at position (" << level << ", "
                  << level_index << ") was invalid." << std::endl;
#endif
        return false;
      }
    }
  }
  return true;
}

template<typename T>
void RandomAccessHeap<T>::Reset(const std::vector<T>& values) {
  values_ = values;
  num_comparison_levels_ = CeilLog2(values_.size());

  // Initialize all of the values as valid.
  valid_values_.clear();
  valid_values_.resize(values_.size(), true);

  if (values_.empty()) {
    return;
  }

  // Compute the comparison metadata from the leaves up (in linear time).
  comparison_tree_.resize(values_.size() - 1);
  for (Int level = num_comparison_levels_ - 1; level >= 0; --level) {
    const Int level_size = LevelSize(level);
    for (Int level_index = 0; level_index < level_size; ++level_index) {
      UpdateComparisonUsingChildren(level, level_index);
    }
  }
}

template<typename T>
std::pair<Int, T> RandomAccessHeap<T>::MinimalEntry() const {
  std::pair<Int, T> entry;
  if (values_.empty()) {
    entry.first = -1;
    entry.second = 0;
    return entry;
  }
  if (values_.size() == 1) {
    entry.first = 0;
    entry.second = values_[entry.first];
    return entry;
  }

  entry.first = comparison_tree_[0];
  entry.second = values_[entry.first];
  return entry;
}

template<typename T>
std::pair<Int, Int> RandomAccessHeap<T>::ParentLevelAndIndex(Int index) const {
  std::pair<Int, Int> tree_pos;
  const Int last_level_size = LevelSize(num_comparison_levels_ - 1);
  if (index < 2 * last_level_size) {
    tree_pos.first = num_comparison_levels_ - 1;
    tree_pos.second = index / 2;
  } else {
    tree_pos.first = num_comparison_levels_ - 2;
    tree_pos.second = (index - last_level_size) / 2;
  }
  return tree_pos;
}

template<typename T>
void RandomAccessHeap<T>::PropagateComparisons(Int index) {
  if (values_.size() <= 1) {
    return;
  }
  std::pair<Int, Int> tree_pos = ParentLevelAndIndex(index);
  UpdateComparisonUsingChildren(tree_pos.first, tree_pos.second);
  while (tree_pos.first != 0) {
    --tree_pos.first;
    tree_pos.second /= 2;
    UpdateComparisonUsingChildren(tree_pos.first, tree_pos.second);
  }
}

template<typename T>
const std::vector<bool>& RandomAccessHeap<T>::ValidValues() const {
  return valid_values_;
}

template<typename T>
bool RandomAccessHeap<T>::ValidValue(Int index) const {
  return valid_values_[index];
}

template<typename T>
const std::vector<T>& RandomAccessHeap<T>::Values() const {
  return values_;
}

template<typename T>
const T& RandomAccessHeap<T>::Value(Int index) const {
  return values_[index];
}

template<typename T>
void RandomAccessHeap<T>::SetValue(Int index, const T& value) {
  if (values_[index] == value) {
    return;
  }
  values_[index] = value;
  if (valid_values_[index]) {
    PropagateComparisons(index);
  }
}

template<typename T>
void RandomAccessHeap<T>::UpdateValue(Int index, const T& value) {
  if (value == T(0)) {
    return;
  }
  values_[index] += value;
  if (valid_values_[index]) {
    PropagateComparisons(index);
  }
}

template<typename T>
void RandomAccessHeap<T>::DisableIndex(Int index) {
  if (!valid_values_[index]) {
    return;
  }
  valid_values_[index] = false;
  PropagateComparisons(index);
}

template<typename T>
Int RandomAccessHeap<T>::PowerOfTwo(UInt exponent) {
  return UInt(1) << exponent;
}

template<typename T>
UInt RandomAccessHeap<T>::CeilLog2(UInt n) {
  UInt ceil_log2 = 0;
  while (n) {
    ++ceil_log2;
    n >>= UInt(1);
  }
  return ceil_log2;
}

} // namespace quotient

#endif // ifndef QUOTIENT_RANDOM_ACCESS_HEAP_H_
