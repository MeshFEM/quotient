/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_RANDOM_ACCESS_HEAP_H_
#define QUOTIENT_RANDOM_ACCESS_HEAP_H_

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
  void UpdateValue(Int index, const T& value_update);

  // Allocates space so that up to 'max_value_changes' calls to
  // 'QueueValueUpdate'/'QueueValueAssignment'/'QueueIndexDisablement' can be
  // performed without another memory allocation.
  void ReserveValueChanges(Int max_value_updates);

  // Modifies the given index to have the specified value but delays updating
  // the comparison tree until the next 'FlushValueChangeQueue' call.
  void QueueValueAssignment(Int index, const T& value);

  // Applies a specified update to the given index to the current update list
  // but delays updating the comparison tree until the next
  // 'FlushValueChangeQueue' call.
  void QueueValueUpdate(Int index, const T& value_update);

  // Marks the given index for disablement but delays the updates to the
  // comparison tree until the next 'FlushValueChangeQueue' call.
  void QueueIndexDisablement(Int index);

  // Updates the comparison tree to reflect the queued value changes.
  void FlushValueChangeQueue();

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

   // The list of indices whose values were updated but whose ancestors in the
   // comparison tree have not yet been updated.
   std::vector<Int> indices_to_update_;

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
};

} // namespace quotient

#include "quotient/random_access_heap-impl.hpp"

#endif // ifndef QUOTIENT_RANDOM_ACCESS_HEAP_H_
