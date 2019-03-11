/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#define CATCH_CONFIG_MAIN
#include <vector>
#include "catch2/catch.hpp"

#include "quotient.hpp"

using quotient::Int;

// Tests for the equality (ignoring set ordering) of a tuple of sets and their
// concatenation in a vector.
template <typename T>
bool SetTuplesAreEqual(const quotient::Buffer<quotient::Buffer<T>>& set_tuple,
                       const quotient::Buffer<T>& concatenation) {
  Int offset = 0;
  for (const quotient::Buffer<T>& set : set_tuple) {
    // Get a sorted copy of this set.
    quotient::Buffer<T> set_copy(set);
    std::sort(set_copy.begin(), set_copy.end());

    // Get a sorted copy of the corresponding portion of the set concatenations.
    const Int set_size = set.Size();
    quotient::Buffer<T> concat_copy(concatenation.begin() + offset,
                                    concatenation.begin() + offset + set_size);
    std::sort(concat_copy.begin(), concat_copy.end());

    // Compare the two sets.
    if (set_copy != concat_copy) {
      return false;
    }

    offset += set_size;
  }
  return true;
}

/* A reproduction of Figs. 1 and 2 from [ADD-96].

   The Amestoy external degree bound produces an overestimate for index
   5 (6 with 1-based indexing) just before the third pivot is selected.
   This helps lead to a deviation between Figs. 1-2 and the AMD results.
*/
TEST_CASE("ADD-96 Figures 1-2", "[ADD-96 Figs 1-2]") {
  quotient::CoordinateGraph graph;
  graph.Resize(10);

  graph.AddEdge(0, 3);
  graph.AddEdge(0, 5);
  graph.AddEdge(3, 0);
  graph.AddEdge(5, 0);

  graph.AddEdge(1, 4);
  graph.AddEdge(1, 5);
  graph.AddEdge(1, 8);
  graph.AddEdge(4, 1);
  graph.AddEdge(5, 1);
  graph.AddEdge(8, 1);

  graph.AddEdge(2, 4);
  graph.AddEdge(2, 5);
  graph.AddEdge(2, 6);
  graph.AddEdge(4, 2);
  graph.AddEdge(5, 2);
  graph.AddEdge(6, 2);

  graph.AddEdge(3, 6);
  graph.AddEdge(3, 7);
  graph.AddEdge(6, 3);
  graph.AddEdge(7, 3);

  graph.AddEdge(4, 6);
  graph.AddEdge(4, 8);
  graph.AddEdge(6, 4);
  graph.AddEdge(8, 4);

  graph.AddEdge(6, 7);
  graph.AddEdge(6, 8);
  graph.AddEdge(6, 9);
  graph.AddEdge(7, 6);
  graph.AddEdge(8, 6);
  graph.AddEdge(9, 6);

  graph.AddEdge(7, 8);
  graph.AddEdge(7, 9);
  graph.AddEdge(8, 7);
  graph.AddEdge(9, 7);

  graph.AddEdge(8, 9);
  graph.AddEdge(9, 8);

  quotient::MinimumDegreeControl control;
  control.degree_type = quotient::kExactDegree;
  control.force_minimal_pivot_indices = true;
  quotient::QuotientGraph quotient_graph(graph, control);
  quotient::MinimumDegreeResult analysis =
      quotient::MinimumDegree(&quotient_graph);

  /* Because of the ordering of the hash bucket, we will prefer the last member
     as the key. But there are several equally-valid solutions (the
     only nontrivial supervariable should be {6, 7, 8}, but the principal
     member can vary).

     The expected supernodes are:

      {0}, {1}, {2}, {3}, {4}, {5}, {}, {}, {6, 7, 8}, {9}.
  */
  const std::vector<Int> kExpectedEliminationOrder{
      0, 1, 2, 3, 4, 5, 8, 9,
  };
  const std::vector<Int> elimination_order = quotient_graph.EliminationOrder();
  REQUIRE(elimination_order == kExpectedEliminationOrder);

  /* This structure is defined directly (modulo translation from 1-based to
     0-based indexing) from the bottom-right of Fig. 2 of [ADD-96]:

      {3, 5}, {4, 5, 8}, {4, 5, 6}, {5, 6, 7}, {5, 6, 8}, {6, 7, 8}, {9}, {}.

     The parents array is thus:

      3, 4, 4, 5, 5, 6, 9, ROOT.

     The resulting assembly tree is:

                    9
                    |
                {6, 7, 8}
                    |
                    5
                   / \
                 3     4
                 |    / \
                 0   1   2

     This implies a lexicographic postordering of:

       {0}, {3}, {1}, {2}, {4}, {5}, {6, 7, 8}, {9}.
  */
  const quotient::Buffer<quotient::Buffer<Int>> kExpectedInversePermutation{
      {0}, {3}, {1}, {2}, {4}, {5}, {6, 7, 8}, {9},
  };
  quotient::Buffer<Int> inverse_permutation;
  quotient_graph.ComputePostorder(&inverse_permutation);
  REQUIRE(SetTuplesAreEqual(kExpectedInversePermutation, inverse_permutation));

  /* The inverse of the lexicographic postordering is:

      {0}, {2}, {3}, {1}, {4}, {5}, {6, 7, 8}, {9}.

    The corresponding permuted supernodal assembly tree is:

                     7
                     |
                     6
                     |
                     5
                    / \
                  1     4
                  |    / \
                  0   2   3

    so the permuted supernodal assembly parents are:

      1, 5, 4, 4, 5, 6, 7, ROOT.
  */
  const quotient::Buffer<quotient::Buffer<Int>> kExpectedPermutation{
      {0}, {2}, {3}, {1}, {4}, {5}, {6, 7, 8}, {9},
  };
  quotient::Buffer<Int> permutation;
  quotient::InvertPermutation(inverse_permutation, &permutation);
  REQUIRE(SetTuplesAreEqual(kExpectedPermutation, permutation));

  const quotient::Buffer<Int> kExpectedPermutedSupernodeSizes{
      1, 1, 1, 1, 1, 1, 3, 1,
  };
  quotient::Buffer<Int> supernode_sizes;
  quotient_graph.PermutedSupernodeSizes(inverse_permutation, &supernode_sizes);
  REQUIRE(kExpectedPermutedSupernodeSizes == supernode_sizes);

  const quotient::Buffer<Int> kExpectedPermutedAssemblyParents{
      1, 5, 4, 4, 5, 6, 7, -1,
  };
  quotient::Buffer<Int> member_to_supernode;
  quotient_graph.PermutedMemberToSupernode(inverse_permutation,
                                           &member_to_supernode);
  quotient::Buffer<Int> parents;
  quotient_graph.PermutedAssemblyParents(permutation, member_to_supernode,
                                         &parents);
  REQUIRE(kExpectedPermutedAssemblyParents == parents);

  const Int kExpectedNumAggressiveAbsorptions = 0;
  REQUIRE(analysis.num_aggressive_absorptions ==
          kExpectedNumAggressiveAbsorptions);
}

// Please see the beginning of Section 5 of [ADD-96].
TEST_CASE("ADD-96 Aggressive Absorbtion", "[ADD-96-Agg-Aborb]") {
  quotient::CoordinateGraph graph;
  graph.Resize(4);

  graph.AddEdge(0, 2);
  graph.AddEdge(0, 3);
  graph.AddEdge(2, 0);
  graph.AddEdge(3, 0);

  graph.AddEdge(1, 2);
  graph.AddEdge(1, 3);
  graph.AddEdge(2, 1);
  graph.AddEdge(3, 1);

  quotient::MinimumDegreeControl control;
  control.degree_type = quotient::kExactDegree;
  control.allow_supernodes = false;
  control.force_minimal_pivot_indices = true;
  control.aggressive_absorption = true;
  quotient::QuotientGraph quotient_graph(graph, control);
  quotient::MinimumDegreeResult analysis =
      quotient::MinimumDegree(&quotient_graph);

  const std::vector<Int> kExpectedEliminationOrder{
      0, 1, 2, 3,
  };
  const std::vector<Int> elimination_order = quotient_graph.EliminationOrder();
  REQUIRE(elimination_order == kExpectedEliminationOrder);

  /* This structure is defined directly (modulo translation from 1-based to
     0-based indexing) from the bottom-right of Fig. 2 of [ADD-96]. The
     elimination structures are:

      {2, 3}, {2, 3}, {3}, {}.

     The resulting parents array is:

      2, 2, 3, ROOT.

     Thus, the assembly tree would be:

        3
        |
        2
       / \
      0   1

     but aggressive absorption causes 1 to be absorbed into 0. The result is an
     assembly tree

        3
        |
        2
        |
        1
        |
        0

     This implies a lexicographic postordering of:

     0, 1, 2, 3.
  */
  const quotient::Buffer<quotient::Buffer<Int>> kExpectedInversePermutation{
      {0}, {1}, {2}, {3},
  };
  quotient::Buffer<Int> inverse_permutation;
  quotient_graph.ComputePostorder(&inverse_permutation);
  REQUIRE(SetTuplesAreEqual(kExpectedInversePermutation, inverse_permutation));

  /* The inverse of the lexicographic postordering is:

      {0}, {1}, {2}, {3}.

    The corresponding permuted supernodal assembly tree is:

        3
        |
        2
        |
        1
        |
        0

    so the permuted supernodal assembly parents are:

      1, 2, 3, ROOT.
  */
  const quotient::Buffer<quotient::Buffer<Int>> kExpectedPermutation{
      {0}, {1}, {2}, {3},
  };
  quotient::Buffer<Int> permutation;
  quotient::InvertPermutation(inverse_permutation, &permutation);
  REQUIRE(SetTuplesAreEqual(kExpectedPermutation, permutation));

  const quotient::Buffer<Int> kExpectedPermutedSupernodeSizes{
      1, 1, 1, 1,
  };
  quotient::Buffer<Int> supernode_sizes;
  quotient_graph.PermutedSupernodeSizes(inverse_permutation, &supernode_sizes);
  REQUIRE(kExpectedPermutedSupernodeSizes == supernode_sizes);

  const quotient::Buffer<Int> kExpectedPermutedAssemblyParents{
      1, 2, 3, -1,
  };
  quotient::Buffer<Int> member_to_supernode;
  quotient_graph.PermutedMemberToSupernode(inverse_permutation,
                                           &member_to_supernode);
  quotient::Buffer<Int> parents;
  quotient_graph.PermutedAssemblyParents(permutation, member_to_supernode,
                                         &parents);
  REQUIRE(kExpectedPermutedAssemblyParents == parents);

  // [ADD-96] discusses the aggressive absorption, 0 into 1.
  const Int kExpectedNumAggressiveAbsorptions = 1;
  REQUIRE(analysis.num_aggressive_absorptions ==
          kExpectedNumAggressiveAbsorptions);
}

// Please see the beginning of Section 5 of [ADD-96].
TEST_CASE("ADD-96 No Aggressive Absorbtion", "[ADD-96-No-Agg-Aborb]") {
  quotient::CoordinateGraph graph;
  graph.Resize(4);

  graph.AddEdge(0, 2);
  graph.AddEdge(0, 3);
  graph.AddEdge(2, 0);
  graph.AddEdge(3, 0);

  graph.AddEdge(1, 2);
  graph.AddEdge(1, 3);
  graph.AddEdge(2, 1);
  graph.AddEdge(3, 1);

  quotient::MinimumDegreeControl control;
  control.degree_type = quotient::kExactDegree;
  control.allow_supernodes = false;
  control.force_minimal_pivot_indices = true;
  control.aggressive_absorption = false;
  quotient::QuotientGraph quotient_graph(graph, control);
  quotient::MinimumDegreeResult analysis =
      quotient::MinimumDegree(&quotient_graph);

  const std::vector<Int> kExpectedEliminationOrder{
      0, 1, 2, 3,
  };
  const std::vector<Int> elimination_order = quotient_graph.EliminationOrder();
  REQUIRE(elimination_order == kExpectedEliminationOrder);

  /* This structure is defined directly (modulo translation from 1-based to
     0-based indexing) from the bottom-right of Fig. 2 of [ADD-96]. The
     elimination structures are:

      {2, 3}, {2, 3}, {3}, {}.

     The resulting parents array is:

      2, 2, 3, ROOT.

     Thus, the assembly tree is:

        3
        |
        2
       / \
      0   1

     This implies a lexicographic postordering of:

     0, 1, 2, 3.
  */
  const quotient::Buffer<quotient::Buffer<Int>> kExpectedInversePermutation{
      {0}, {1}, {2}, {3},
  };
  quotient::Buffer<Int> inverse_permutation;
  quotient_graph.ComputePostorder(&inverse_permutation);
  REQUIRE(SetTuplesAreEqual(kExpectedInversePermutation, inverse_permutation));

  /* The inverse of the lexicographic postordering is:

      {0}, {1}, {2}, {3}.

    The corresponding permuted supernodal assembly tree is:

        3
        |
        2
       / \
      0   1

    so the permuted supernodal assembly parents are:

      2, 2, 3, ROOT.
  */
  const quotient::Buffer<quotient::Buffer<Int>> kExpectedPermutation{
      {0}, {1}, {2}, {3},
  };
  quotient::Buffer<Int> permutation;
  quotient::InvertPermutation(inverse_permutation, &permutation);
  REQUIRE(SetTuplesAreEqual(kExpectedPermutation, permutation));

  const quotient::Buffer<Int> kExpectedPermutedSupernodeSizes{
      1, 1, 1, 1,
  };
  quotient::Buffer<Int> supernode_sizes;
  quotient_graph.PermutedSupernodeSizes(inverse_permutation, &supernode_sizes);
  REQUIRE(kExpectedPermutedSupernodeSizes == supernode_sizes);

  const quotient::Buffer<Int> kExpectedPermutedAssemblyParents{
      2, 2, 3, -1,
  };
  quotient::Buffer<Int> member_to_supernode;
  quotient_graph.PermutedMemberToSupernode(inverse_permutation,
                                           &member_to_supernode);
  quotient::Buffer<Int> parents;
  quotient_graph.PermutedAssemblyParents(permutation, member_to_supernode,
                                         &parents);
  REQUIRE(kExpectedPermutedAssemblyParents == parents);

  const Int kExpectedNumAggressiveAbsorptions = 0;
  REQUIRE(analysis.num_aggressive_absorptions ==
          kExpectedNumAggressiveAbsorptions);
}
