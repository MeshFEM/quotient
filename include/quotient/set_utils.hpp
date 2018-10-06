/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_SET_UTILS_H_
#define QUOTIENT_SET_UTILS_H_

#include <iostream>
#include <vector>

#include "quotient/config.hpp"

namespace quotient {

// Returns the number of entries in {vec0} \cap {vec1}, where both vectors
// are assumed sorted and unique.
template<typename T>
Int SizeOfIntersection(
    const std::vector<T>& vec0, const std::vector<T>& vec1) {
  Int num_intersections = 0;
  auto vec0_iter = vec0.cbegin(); 
  auto vec1_iter = vec1.cbegin();
  while (vec0_iter != vec0.cend() && vec1_iter != vec1.cend()) {
    if (*vec0_iter < *vec1_iter) {
      vec0_iter = std::lower_bound(vec0_iter, vec0.cend(), *vec1_iter);
    } else if (*vec0_iter > *vec1_iter) {
      vec1_iter = std::lower_bound(vec1_iter, vec1.cend(), *vec0_iter);
    } else {
      ++num_intersections;
      ++vec0_iter;
      ++vec1_iter;
    }
  }
  return num_intersections;
}

// Returns the number of entries in {vec} \ {blacklist}, where both vectors are
// assumed to be sorted and unique.
template<typename T>
Int SizeOfDifference(
    const std::vector<T>& vec, const std::vector<T>& blacklist) {
  return vec.size() - SizeOfIntersection(vec, blacklist);
}

// Returns the number of entries in {vec0} \cup {vec1}, where both vectors are
// assumed to be sorted and unique.
template<typename T>
Int SizeOfUnion(const std::vector<T>& vec0, const std::vector<T>& vec1) {
  return vec0.size() + vec1.size() - SizeOfIntersection(vec0, vec1);
}

// Returns the number of entries in ({vec0} \cap {vec1}) \ {blacklist}, where
// all vectors are assumed sorted and unique.
template<typename T>
Int SizeOfBlacklistedIntersection(
    const std::vector<T>& vec0,
    const std::vector<T>& vec1,
    const std::vector<T>& blacklist) {
  auto vec0_iter = vec0.cbegin(); 
  auto blacklist_iter0 = blacklist.cbegin();

  auto vec1_iter = vec1.cbegin();
  auto blacklist_iter1 = blacklist.cbegin();

  Int num_blacklisted_intersections = 0;
  while (vec0_iter != vec0.cend() && vec1_iter != vec1.cend()) {
    // Skip this entry of vec0 if it is blacklisted.
    blacklist_iter0 =
      std::lower_bound(blacklist_iter0, blacklist.cend(), *vec0_iter);
    if (blacklist_iter0 != blacklist.cend() && *vec0_iter == *blacklist_iter0) {
      ++vec0_iter; 
      ++blacklist_iter0;
      continue;
    }

    // Skip this entry of vec1 if it is blacklisted.
    blacklist_iter1 =
      std::lower_bound(blacklist_iter1, blacklist.cend(), *vec1_iter);
    if (blacklist_iter1 != blacklist.cend() && *vec1_iter == *blacklist_iter1) {
      ++vec1_iter;
      ++blacklist_iter1;
      continue;
    }

    if (*vec0_iter < *vec1_iter) {
      vec0_iter = std::lower_bound(vec0_iter, vec0.cend(), *vec1_iter);
    } else if (*vec0_iter > *vec1_iter) {
      vec1_iter = std::lower_bound(vec1_iter, vec1.cend(), *vec0_iter);
    } else {
      ++num_blacklisted_intersections;
      ++vec0_iter;
      ++vec1_iter;
    }
  }
  return num_blacklisted_intersections;
}

// Returns the number of entries in ({vec0} \cup {vec1}) \ {blacklist}, where
// all vectors are assumed sorted and unique.
template<typename T>
Int SizeOfBlacklistedUnion(
    const std::vector<T>& vec0,
    const std::vector<T>& vec1,
    const std::vector<T>& blacklist) {
  const Int num_vec0_blacklisted = SizeOfIntersection(vec0, blacklist);
  const Int num_vec1_blacklisted = SizeOfIntersection(vec1, blacklist);
  const Int num_blacklisted_intersections = SizeOfBlacklistedIntersection(
      vec0, vec1, blacklist);
  return (vec0.size() - num_vec0_blacklisted) +
      (vec1.size() - num_vec1_blacklisted) - num_blacklisted_intersections;
}

// Fills 'filtered_vec' with the sorted serialization of {vec} \ {blacklist},
// where 'vec' and 'blacklist' are both sorted, unique vectors.
template<typename T>
void FilterSet(
    const std::vector<T>& vec,
    const std::vector<T>& blacklist,
    std::vector<T>* filtered_vec) {
  const Int filtered_size = SizeOfDifference(vec, blacklist);
  filtered_vec->resize(0);
  filtered_vec->resize(filtered_size);
#ifdef QUOTIENT_DEBUG
  auto iter = std::set_difference(
      vec.begin(), vec.end(),
      blacklist.begin(), blacklist.end(),
      filtered_vec->begin());
  if (filtered_size != std::distance(filtered_vec->begin(), iter)) {
    std::cerr << "Filtered sizes did not match in FilterSet." << std::endl;
  }
#else
  std::set_difference(
      vec.begin(), vec.end(),
      blacklist.begin(), blacklist.end(),
      filtered_vec->begin());
#endif
}

// Fills 'sorted_union' with the sorted serialization of {vec0} \cup {vec1},
// where 'vec0' and 'vec1' are sorted, unique vectors.
template<typename T>
void MergeSets(
    const std::vector<T>& vec0,
    const std::vector<T>& vec1, 
    std::vector<T>* sorted_union) {
  const Int union_size = SizeOfUnion(vec0, vec1);
  sorted_union->resize(0);
  sorted_union->resize(union_size);
#ifdef QUOTIENT_DEBUG
  auto iter = std::set_union(
      vec0.begin(), vec0.end(),
      vec1.begin(), vec1.end(),
      sorted_union->begin());
  if (union_size != std::distance(sorted_union->begin(), iter)) {
    std::cerr << "Union sizes did not match in MergeSets." << std::endl;
  }
#else
  std::set_union(
      vec0.begin(), vec0.end(),
      vec1.begin(), vec1.end(),
      sorted_union->begin());
#endif
}

// Inserts 'value' into the sorted, unique vector, 'vec', so that the result
// is again sorted and unique.
template<typename T>
void InsertNewEntryIntoSet(const T& value, std::vector<T>* vec) {
  auto iter = std::lower_bound(vec->begin(), vec->end(), value);
#ifdef QUOTIENT_DEBUG
  if (iter != vec->end() && *iter == value) {
    std::cerr << value << " was already in set in InsertEntryIntoSet."
              << std::endl;
  }
#endif
  vec->insert(iter, value);
}

// Inserts 'value' into the sorted, unique vector, 'vec', so that the result
// is again sorted and unique.
template<typename T>
void InsertEntryIntoSet(const T& value, std::vector<T>* vec) {
  auto iter = std::lower_bound(vec->begin(), vec->end(), value);
  if (iter == vec->end() || *iter != value) {
    vec->insert(iter, value);
  }
}

} // namespace quotient

#endif // ifndef QUOTIENT_SET_UTILS_H_
