/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_SET_UTILS_H_
#define QUOTIENT_SET_UTILS_H_

#include <vector>

#include "quotient/config.hpp"

namespace quotient {

// Returns the number of entries in {vec0} \cap {vec1}, where both vectors
// are assumed sorted and unique.
//
// TODO(Jack Poulson): Benchmark several different variants of this routine.
template<typename T>
Int SizeOfIntersection(
    const std::vector<T>& vec0, const std::vector<T>& vec1);

// Returns the number of entries in {vec} \ {blacklist}, where both vectors are
// assumed to be sorted and unique.
template<typename T>
Int SizeOfDifference(
    const std::vector<T>& vec, const std::vector<T>& blacklist);

// Returns the number of entries in {vec0} \cup {vec1}, where both vectors are
// assumed to be sorted and unique.
template<typename T>
Int SizeOfUnion(const std::vector<T>& vec0, const std::vector<T>& vec1);

// Returns the number of entries in ({vec0} \cap {vec1}) \ {blacklist}, where
// all vectors are assumed sorted and unique.
template<typename T>
Int SizeOfBlacklistedIntersection(
    const std::vector<T>& vec0,
    const std::vector<T>& vec1,
    const std::vector<T>& blacklist);

// Returns the number of entries in ({vec0} \cup {vec1}) \ {blacklist}, where
// all vectors are assumed sorted and unique.
template<typename T>
Int SizeOfBlacklistedUnion(
    const std::vector<T>& vec0,
    const std::vector<T>& vec1,
    const std::vector<T>& blacklist);

// Fills 'filtered_vec' with the sorted serialization of {vec} \ {blacklist},
// where 'vec' and 'blacklist' are both sorted, unique vectors.
template<typename T>
void FilterSet(
    const std::vector<T>& vec,
    const std::vector<T>& blacklist,
    std::vector<T>* filtered_vec);

// Overwrites 'vec' with {vec} \ {blacklist}.
template<typename T>
void FilterSetInPlace(const std::vector<T>& blacklist, std::vector<T>* vec);

// Overwrites 'vec' with {vec} \ {blacklist0 \cup blacklist1}.
template<typename T>
void DoubleFilterSetInPlace(
    const std::vector<T>& blacklist0,
    const std::vector<T>& blacklist1,
    std::vector<T>* vec);

// Fills 'sorted_union' with the sorted serialization of {vec0} \cup {vec1},
// where 'vec0' and 'vec1' are sorted, unique vectors.
template<typename T>
void MergeSets(
    const std::vector<T>& vec0,
    const std::vector<T>& vec1, 
    std::vector<T>* sorted_union);

// Inserts 'value' into the sorted, unique vector, 'vec', so that the result
// is again sorted and unique.
template<typename T>
void InsertNewEntryIntoSet(const T& value, std::vector<T>* vec);

// Inserts 'value' into the sorted, unique vector, 'vec', so that the result
// is again sorted and unique.
template<typename T>
void InsertEntryIntoSet(const T& value, std::vector<T>* vec);

} // namespace quotient

#include "quotient/set_utils-impl.hpp"

#endif // ifndef QUOTIENT_SET_UTILS_H_
