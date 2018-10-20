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

TEST_CASE("Simple example", "[simple]") {
  using quotient::Int;

  const std::vector<double> orig_values{
      3., 1., 4., 1., 5., 9., 2., 6., 5., 3., 5., 8., 9., 7., 9., 3., 2., 3.,
      8., 4., 6., 2., 6., 4., 3., 3., 8., 3., 2., 7., 9., 5., 0., 2., 8., 8.,
  };

  quotient::RandomAccessHeap<double> heap; 
  heap.Reset(orig_values);

  // The heap now contains:
  //  3., 1., 4., 1., 5., 9., 2., 6., 5., 3., 5., 8., 9., 7., 9., 3., 2., 3.,
  //  8., 4., 6., 2., 6., 4., 3., 3., 8., 3., 2., 7., 9., 5., 0., 2., 8., 8.,
  // with the unique minimum being 0 in position 32.

  std::pair<Int, double> minimal_entry = heap.MinimalEntry();

  REQUIRE(minimal_entry.first == 32);
  REQUIRE(minimal_entry.second == 0.);

  heap.SetValue(7, 3.);
  minimal_entry = heap.MinimalEntry();

  // The heap now contains:
  //  3., 1., 4., 1., 5., 9., 2., 3., 5., 3., 5., 8., 9., 7., 9., 3., 2., 3.,
  //  8., 4., 6., 2., 6., 4., 3., 3., 8., 3., 2., 7., 9., 5., 0., 2., 8., 8.,
  // with the unique minimum being 0 in position 32.

  REQUIRE(minimal_entry.first == 32);
  REQUIRE(minimal_entry.second == 0.);

  heap.SetValue(6, -1.);
  minimal_entry = heap.MinimalEntry();

  // The heap now contains:
  //  3., 1., 4., 1., 5., 9., -1., 3., 5., 3., 5., 8., 9., 7., 9., 3., 2., 3.,
  //  8., 4., 6., 2., 6., 4.,  3., 3., 8., 3., 2., 7., 9., 5., 0., 2., 8., 8.,
  // with the unique minimum being -1 in position 6.

  REQUIRE(minimal_entry.first == 6);
  REQUIRE(minimal_entry.second == -1.);

  heap.DisableIndex(6);
  minimal_entry = heap.MinimalEntry();

  // The heap now contains:
  //  3., 1., 4., 1., 5., 9.,   , 3., 5., 3., 5., 8., 9., 7., 9., 3., 2., 3.,
  //  8., 4., 6., 2., 6., 4., 3., 3., 8., 3., 2., 7., 9., 5., 0., 2., 8., 8.,
  // with the unique minimum being 0 in position 32.

  REQUIRE(minimal_entry.first == 32);
  REQUIRE(minimal_entry.second == 0.);

  for (Int i = 10; i < 20; ++i) {
    heap.DisableIndex(i);
  }
  minimal_entry = heap.MinimalEntry();

  // The heap now contains:
  //  3., 1., 4., 1., 5., 9.,   , 3., 5., 3.,   ,   ,   ,   ,   ,   ,   ,   ,
  //    ,   , 6., 2., 6., 4., 3., 3., 8., 3., 2., 7., 9., 5., 0., 2., 8., 8.,
  // with the unique minimum being 0 in position 32.
  
  REQUIRE(minimal_entry.first == 32);
  REQUIRE(minimal_entry.second == 0.);

  heap.DisableIndex(33);
  minimal_entry = heap.MinimalEntry();

  // The heap now contains:
  //  3., 1., 4., 1., 5., 9.,   , 3., 5., 3.,   ,   ,   ,   ,   ,   ,   ,   ,
  //    ,   , 6., 2., 6., 4., 3., 3., 8., 3., 2., 7., 9., 5., 0.,   , 8., 8.,
  // with the unique minimum being 0 in position 32.

  REQUIRE(minimal_entry.first == 32);
  REQUIRE(minimal_entry.second == 0.);

  heap.DisableIndex(32);
  minimal_entry = heap.MinimalEntry();

  // The heap now contains:
  //  3., 1., 4., 1., 5., 9.,   , 3., 5., 3.,   ,   ,   ,   ,   ,   ,   ,   ,
  //    ,   , 6., 2., 6., 4., 3., 3., 8., 3., 2., 7., 9., 5.,   ,   , 8., 8.,
  // with minima of 1 appearing in positions 1 and 3.

  REQUIRE((minimal_entry.first == 1 || minimal_entry.first == 3));
  REQUIRE(minimal_entry.second == 1.);

  heap.DisableIndex(30);
  heap.DisableIndex(28);
  heap.UpdateValue(2, -6);
  heap.DisableIndex(4);
  heap.UpdateValue(3, -2);
  minimal_entry = heap.MinimalEntry();

  // The heap now contains:
  //  3., 1., -2., -1.,   , 9.,   , 3., 5., 3.,   ,   ,   ,   ,   ,   ,   ,   ,
  //    ,   ,  6.,  2., 6., 4., 3., 3., 8., 3.,   , 7.,   , 5.,   ,   , 8., 8.,
  // with minimum of -2 in position 2.

  REQUIRE(minimal_entry.first == 2);
  REQUIRE(minimal_entry.second == -2.);
}

TEST_CASE("Disable last tree index", "[disable-last]") {
  using quotient::Int;

  const std::vector<double> orig_values{
      3., 1., 4., 1., 5., 9., 2., 6., 5., 3., 5., 8., 9., 7., 9., 3., 2.,  3.,
      8., 4., 6., 2., 6., 4., 3., 3., 8., 3., 2., 7., 9., 5., 0., 2., 8., -1.,
  };

  quotient::RandomAccessHeap<double> heap; 
  heap.Reset(orig_values);

  // The heap now contains:
  //  3., 1., 4., 1., 5., 9., 2., 6., 5., 3., 5., 8., 9., 7., 9., 3., 2.,  3.,
  //  8., 4., 6., 2., 6., 4., 3., 3., 8., 3., 2., 7., 9., 5., 0., 2., 8., -1.,
  // with the unique minimum being -1 in position 35.

  std::pair<Int, double> minimal_entry = heap.MinimalEntry();

  REQUIRE(minimal_entry.first == 35);
  REQUIRE(minimal_entry.second == -1.);

  heap.DisableIndex(35);
  minimal_entry = heap.MinimalEntry();

  // The heap now contains:
  //  3., 1., 4., 1., 5., 9., 2., 6., 5., 3., 5., 8., 9., 7., 9., 3., 2.,  3.,
  //  8., 4., 6., 2., 6., 4., 3., 3., 8., 3., 2., 7., 9., 5., 0., 2., 8.,    ,
  // with the unique minimum being 0 in position 32.

  REQUIRE(minimal_entry.first == 32);
  REQUIRE(minimal_entry.second == 0.);
}

TEST_CASE("Batch example", "[batch]") {
  using quotient::Int;

  const std::vector<double> orig_values{
      3., 1., 4., 1., 5., 9., 2., 6., 5., 3., 5., 8., 9., 7., 9., 3., 2., 3.,
      8., 4., 6., 2., 6., 4., 3., 3., 8., 3., 2., 7., 9., 5., 0., 2., 8., 8.,
  };

  quotient::RandomAccessHeap<double> heap; 
  heap.Reset(orig_values);

  // The heap now contains:
  //  3., 1., 4., 1., 5., 9., 2., 6., 5., 3., 5., 8., 9., 7., 9., 3., 2., 3.,
  //  8., 4., 6., 2., 6., 4., 3., 3., 8., 3., 2., 7., 9., 5., 0., 2., 8., 8.,
  // with the unique minimum being 0 in position 32.

  std::pair<Int, double> minimal_entry = heap.MinimalEntry();

  REQUIRE(minimal_entry.first == 32);
  REQUIRE(minimal_entry.second == 0.);

  heap.ReserveValueChanges(20);

  heap.QueueValueAssignment(7, 3.);
  heap.QueueValueAssignment(6, -1.);
  heap.QueueIndexDisablement(6);
  for (Int i = 10; i < 20; ++i) {
    heap.QueueIndexDisablement(i);
  }
  heap.QueueIndexDisablement(33);
  heap.QueueIndexDisablement(32);
  heap.QueueIndexDisablement(30);
  heap.QueueIndexDisablement(28);
  heap.QueueValueUpdate(2, -6);
  heap.QueueIndexDisablement(4);
  heap.QueueValueUpdate(3, -2);

  heap.FlushValueChangeQueue();
  minimal_entry = heap.MinimalEntry();

  // The heap now contains:
  //  3., 1., -2., -1.,   , 9.,   , 3., 5., 3.,   ,   ,   ,   ,   ,   ,   ,   ,
  //    ,   ,  6.,  2., 6., 4., 3., 3., 8., 3.,   , 7.,   , 5.,   ,   , 8., 8.,
  // with minimum of -2 in position 2.

  REQUIRE(minimal_entry.first == 2);
  REQUIRE(minimal_entry.second == -2.);
}
