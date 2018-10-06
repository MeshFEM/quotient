/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#define CATCH_CONFIG_MAIN
#include <vector>
#include "quotient.hpp"
#include "catch.hpp"

TEST_CASE("SizeOfIntersection", "[size-of-intersection]") {
  const std::vector<quotient::Int> vec0{
      0, 10, 12, 13, 18, 25, 30, 72, 73, 74, 75, 81, 95, 115,
  };
  const std::vector<quotient::Int> vec1{
      1, 2, 3, 4, 10, 11, 12, 13, 14, 19, 20, 24, 25, 26, 27, 28, 29, 30, 31,
      70, 71, 72, 74, 75, 76, 77, 81,
  };
  REQUIRE(quotient::SizeOfIntersection(vec0, vec1) == 9);
}

TEST_CASE("SizeOfDifference", "[size-of-difference]") {
  const std::vector<quotient::Int> vec0{
      0, 10, 12, 13, 18, 25, 30, 72, 73, 74, 75, 81, 95, 115,
  };
  const std::vector<quotient::Int> vec1{
      1, 2, 3, 4, 10, 11, 12, 13, 14, 19, 20, 24, 25, 26, 27, 28, 29, 30, 31,
      70, 71, 72, 74, 75, 76, 77, 81,
  };
  REQUIRE(quotient::SizeOfDifference(vec0, vec1) == 5);
}

TEST_CASE("SizeOfUnion", "[size-of-union]") {
  const std::vector<quotient::Int> vec0{
      0, 10, 12, 13, 18, 25, 30, 72, 73, 74, 75, 81, 95, 115,
  };
  const std::vector<quotient::Int> vec1{
      1, 2, 3, 4, 10, 11, 12, 13, 14, 19, 20, 24, 25, 26, 27, 28, 29, 30, 31,
      70, 71, 72, 74, 75, 76, 77, 81,
  };
  REQUIRE(quotient::SizeOfUnion(vec0, vec1) == 32);
}

TEST_CASE("SizeOfBlacklistedIntersection", "[size-of-blacklist-intersect]") {
  const std::vector<quotient::Int> vec0{
      0, 10, 12, 13, 18, 25, 30, 72, 73, 74, 75, 81, 95, 115,
  };
  const std::vector<quotient::Int> vec1{
      1, 2, 3, 4, 10, 11, 12, 13, 14, 19, 20, 24, 25, 26, 27, 28, 29, 30, 31,
      70, 71, 72, 74, 75, 76, 77, 81,
  };
  const std::vector<quotient::Int> blacklist{
      2, 10, 13, 17, 19, 28, 29, 44, 45, 46, 72, 81, 90, 110, 200,
  };
  REQUIRE(quotient::SizeOfBlacklistedIntersection(vec0, vec1, blacklist) == 5);
}

TEST_CASE("SizeOfBlacklistedUnion", "[size-of-blacklist-union]") {
  const std::vector<quotient::Int> vec0{
      0, 10, 12, 13, 18, 25, 30, 72, 73, 74, 75, 81, 95, 115,
  };
  const std::vector<quotient::Int> vec1{
      1, 2, 3, 4, 10, 11, 12, 13, 14, 19, 20, 24, 25, 26, 27, 28, 29, 30, 31,
      70, 71, 72, 74, 75, 76, 77, 81,
  };
  const std::vector<quotient::Int> blacklist{
      2, 10, 13, 17, 19, 28, 29, 44, 45, 46, 72, 81, 90, 110, 200,
  };
  REQUIRE(quotient::SizeOfBlacklistedUnion(vec0, vec1, blacklist) == 24);
}

TEST_CASE("FilterSet", "[filter-set]") {
  const std::vector<quotient::Int> vec{
      0, 10, 12, 13, 18, 25, 30, 72, 73, 74, 75, 81, 95, 115,
  };
  const std::vector<quotient::Int> blacklist{
      2, 10, 13, 17, 19, 28, 29, 44, 45, 46, 72, 81, 90, 110, 200,
  };
  const std::vector<quotient::Int> kExpectedFilteredVec{
      0, 12, 18, 25, 30, 73, 74, 75, 95, 115,
  };

  std::vector<quotient::Int> filtered_vec;
  quotient::FilterSet(vec, blacklist, &filtered_vec);

  REQUIRE(filtered_vec == kExpectedFilteredVec);
}

TEST_CASE("MergeSets", "[merge-sets]") {
  const std::vector<quotient::Int> vec0{
      0, 10, 12, 13, 18, 25, 30, 72, 73, 74, 75, 81, 95, 115,
  };
  const std::vector<quotient::Int> vec1{
      1, 2, 3, 4, 10, 11, 12, 13, 14, 19, 20, 24, 25, 26, 27, 28, 29, 30, 31,
      70, 71, 72, 74, 75, 76, 77, 81,
  };
  const std::vector<quotient::Int> kExpectedUnion{
      0, 1, 2, 3, 4, 10, 11, 12, 13, 14, 18, 19, 20, 24, 25, 26, 27, 28, 29, 30,
      31, 70, 71, 72, 73, 74, 75, 76, 77, 81, 95, 115,
  };

  std::vector<quotient::Int> sorted_union;
  quotient::MergeSets(vec0, vec1, &sorted_union);

  REQUIRE(sorted_union == kExpectedUnion);
}

TEST_CASE("InsertEntryIntoSet", "[insert-into-set]") {
  std::vector<quotient::Int> vec{
      1, 2, 3, 4, 10, 11, 12, 13, 14, 19, 20, 24, 25, 26, 27, 28, 29, 30, 31,
      70, 71, 72, 74, 75, 76, 77, 81,
  };

  quotient::InsertEntryIntoSet(quotient::Int(0), &vec);
  const std::vector<quotient::Int> kExpectedVec0{
      0, 1, 2, 3, 4, 10, 11, 12, 13, 14, 19, 20, 24, 25, 26, 27, 28, 29, 30, 31,
      70, 71, 72, 74, 75, 76, 77, 81,
  };
  REQUIRE(vec == kExpectedVec0);

  quotient::InsertNewEntryIntoSet(quotient::Int(18), &vec);
  const std::vector<quotient::Int> kExpectedVec1{
      0, 1, 2, 3, 4, 10, 11, 12, 13, 14, 18, 19, 20, 24, 25, 26, 27, 28, 29, 30,
      31, 70, 71, 72, 74, 75, 76, 77, 81,
  };
  REQUIRE(vec == kExpectedVec1);

  quotient::InsertEntryIntoSet(quotient::Int(85), &vec);
  const std::vector<quotient::Int> kExpectedVec2{
      0, 1, 2, 3, 4, 10, 11, 12, 13, 14, 18, 19, 20, 24, 25, 26, 27, 28, 29, 30,
      31, 70, 71, 72, 74, 75, 76, 77, 81, 85,
  };
  REQUIRE(vec == kExpectedVec2);

  quotient::InsertEntryIntoSet(quotient::Int(18), &vec);
  REQUIRE(vec == kExpectedVec2);
}
