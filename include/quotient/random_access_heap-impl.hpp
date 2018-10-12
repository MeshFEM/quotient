/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_RANDOM_ACCESS_HEAP_IMPL_H_
#define QUOTIENT_RANDOM_ACCESS_HEAP_IMPL_H_

#include <iostream>
#include <vector>

#include "quotient/config.hpp"
#include "quotient/random_access_heap.hpp"

namespace quotient {

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

#endif // ifndef QUOTIENT_RANDOM_ACCESS_HEAP_IMPL_H_
