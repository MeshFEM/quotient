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
// TODO(Jack Poulson): Add support for a PushBack(const T& value) method and
// a DeleteIndex(Int index) method.
//
// TODO(Jack Poulson): Provide support for overloading operator <.
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
   Int ChildMinimalTreeIndex(
       Int level, Int level_index, bool right_child) const;

   // Updates all relevant ancestral metadata for the leaf node with given
   // index.
   void PropagateComparisons(Int index);

   // Updates the entry of 'comparison_tree_' corresponding to the given level
   // and index within said level using its two children.
   void UpdateComparisonUsingChildren(Int level, Int level_index);

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

   /* An array of length 'num_active_ - 1' that stores the packed indices of
      the minimal (permuted) active index of the minimal descendant.

      For example, if there are four active indices, say, with values:
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
      level is truncated so that there are exactly 'num_active_ - 1' nodes in
      the tree. For example, if there are five active indices, say, with values:
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
        [2^j - 1, ..., max(2^{j + 1} - 2, num_active_ - 1)].
   */
   std::vector<Int> comparison_tree_;

   // The number of levels of comparison metadata.
   Int num_comparison_levels_;
};

template<typename T>
RandomAccessHeap<T>::RandomAccessHeap() { }

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
    return num_active_ - PowerOfTwo(num_comparison_levels_ - 1);
  }
  return PowerOfTwo(level);
}

template<typename T>
Int RandomAccessHeap<T>::ChildMinimalTreeIndex(
    Int level, Int level_index, bool right_child) const {
#ifdef QUOTIENT_DEBUG
  if (level < 0 || level >= num_comparison_levels_) {
    std::cerr << "Requested level size of invalid level." << std::endl;
    return;
  }
  if (level_index < 0 || index >= LevelSize(level)) {
    std::cerr << "Requested invalid index of level." << std::endl;
    return;
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
void RandomAccessHeap<T>::UpdateComparisonUsingChildren(
    Int level, Int level_index) {
  const Int left_tree_index = ChildMinimalTreeIndex(
      level, level_index, false /* right_child */);
  const Int right_tree_index = ChildMinimalTreeIndex(
      level, level_index, true /* right_child */);

  const T& left_value = values_[inverse_perm_[left_tree_index]];
  const T& right_value = values_[inverse_perm_[right_tree_index]];

  const Int new_tree_index =
      left_value <= right_value ? left_tree_index : right_tree_index;

  const Int level_offset = LevelOffset(level);
  comparison_tree_[level_offset + level_index] = new_tree_index;
}

template<typename T>
void RandomAccessHeap<T>::Initialize(const std::vector<T>& values) {
  values_ = values;
  num_active_ = values.size();
  num_comparison_levels_ = CeilLog2(num_active_);
  if (num_active_ == 0) {
    return;
  }

  // Create an identity permutation.
  perm_.resize(values.size());
  std::iota(perm_.begin(), perm_.end(), 0);
  inverse_perm_ = perm_;

  // Compute the comparison metadata from the leaves up (in linear time).
  comparison_tree_.resize(num_active_ - 1);
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
  entry.first = inverse_perm_[comparison_tree_[0]];
  entry.second = values_[entry.first];
  return entry;
}

template<typename T>
std::pair<Int, Int> RandomAccessHeap<T>::ParentLevelAndIndex(
    Int tree_index) const {
  std::pair<Int, Int> tree_pos;
  const Int last_level_size = LevelSize(num_comparison_levels_ - 1);
  if (tree_index < 2 * last_level_size) {
    tree_pos.first = num_comparison_levels_ - 1;
    tree_pos.second = tree_index / 2;
  } else {
    tree_pos.first = num_comparison_levels_ - 2;
    tree_pos.second = (tree_index - last_level_size) / 2;
  }
  return tree_pos;
}

template<typename T>
void RandomAccessHeap<T>::PropagateComparisons(Int index) {
  const Int tree_index = perm_[index];
  std::pair<Int, Int> tree_pos = ParentLevelAndIndex(tree_index);
  UpdateComparisonUsingChildren(tree_pos.first, tree_pos.second);
  while (tree_pos.first != 0) {
    --tree_pos.first;
    tree_pos.second /= 2;
    UpdateComparisonUsingChildren(tree_pos.first, tree_pos.second);
  }
}

template<typename T>
const T& RandomAccessHeap<T>::GetValue(Int index) const {
  return values_[index];
}

template<typename T>
void RandomAccessHeap<T>::SetValue(Int index, const T& value) {
  if (values_[index] == value) {
    return;
  }
  values_[index] = value;
  PropagateComparisons(index);
}

template<typename T>
void RandomAccessHeap<T>::UpdateValue(Int index, const T& value) {
  if (value == T(0)) {
    return;
  }
  values_[index] += value;
  PropagateComparisons(index);
}

template<typename T>
void RandomAccessHeap<T>::DisableIndex(Int index) {
  const Int old_tree_index = perm_[index];
#ifdef QUOTIENT_DEBUG
  if (old_tree_index >= num_active_) {
    std::cerr << "Index was already disabled." << std::endl;
    return;
  }
#endif

  // Apply the transposition (old_tree_index, num_active_ - 1) to the current
  // permutation from the left.
  //
  // More generally, a permutation and its inverse may be quickly updated to
  // represent the original permutation with a swap (a, b) applied from the
  // left by:
  //
  //   1) Storing a' := inverse_perm[a], b' := inverse_perm[b].
  //
  //   2) Updating inverse_perm via swap(inverse_perm[a], inverse_perm[b]).
  //
  //   3) Updating perm via swap(perm[a'], perm[b']).
  //
  // But, since exchanging steps (2) and (1) does not change the effect of the
  // sequence, we may instead execute:
  //
  //   1) Updating inverse_perm via swap(inverse_perm[a], inverse_perm[b]).
  //
  //   2) Updating perm via swap(perm[inverse_perm[a]], perm[inverse_perm[b]]).
  //
  const Int swap_ind0 = old_tree_index;
  const Int swap_ind1 = num_active_ - 1;
  std::swap(inverse_perm_[swap_ind0], inverse_perm_[swap_ind1]);
  std::swap(perm_[inverse_perm_[swap_ind0]], perm_[inverse_perm_[swap_ind1]]);

  // Disable the last index of the tree.
  --num_active_;
  comparison_tree_.pop_back();
  num_comparison_levels_ = CeilLog2(num_active_);

  // Ensure that the comparison metadata is udpated. There are at most two
  // locations where we must manually propagate changes:
  //
  //  1) Tree index 'old_tree_index'.
  //
  //  2) Tree index 'num_active_ - 1'.
  //
  // The former must be skipped if 'old_tree_index == num_active_'.
  //
  if (old_tree_index != num_active_) {
    PropagateComparisons(inverse_perm_[old_tree_index]);
  }
  PropagateComparisons(inverse_perm_[num_active_ - 1]);
}

template<typename T>
inline Int RandomAccessHeap<T>::PowerOfTwo(UInt exponent) {
  return UInt(1) << exponent;
}

template<typename T>
inline UInt RandomAccessHeap<T>::CeilLog2(UInt n) {
  UInt ceil_log2 = 0;
  while (n) {
    ++ceil_log2;
    n >>= UInt(1);
  }
  return ceil_log2;
}

} // namespace quotient

#endif // ifndef QUOTIENT_RANDOM_ACCESS_HEAP_H_
