/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_QUOTIENT_GRAPH_H_
#define QUOTIENT_QUOTIENT_GRAPH_H_

#include <vector>

#include "quotient/config.hpp"
#include "quotient/coordinate_graph.hpp"

namespace quotient {

// A choice of hash function for mapping a variable to an std::size_t.
enum VariableHashType {
  // Due to Ashcraft and described in [ADD-96].
  kAshcraftVariableHash,

  // A custom hash function meant to avoid collisions by removing the
  // modular arithmetic (allowing the full range of std::size_t) and
  // multiplying each index update by its position in the adjaency or
  // element list.
  kBasicVariableHash,
};

// A doubly-linked list for tracking node degrees in a manner which makes it
// easy to find a minimal degree variable.
struct DegreeLists {
  // A lower-bound on the minimum degree of all degree list members.
  Int degree_lower_bound = 0;

  // A list of the external degrees associated with each index.
  std::vector<Int> degrees;

  // The indices of the heads of the degree lists for each degree. If no
  // index exists for degree 'degree', then 'degree_heads[degree]' is -1.
  std::vector<Int> degree_heads;

  // Index 'i' will provide the index of the member that occurs directly after
  // 'i' in the degree list containing index 'i' (if such a member does not
  // exit, the value will be -1).
  std::vector<Int> next_degree_member;

  // Index 'i' will provide the index of the member that occurs directly before
  // 'i' in the degree list containing index 'i' (if such a member does not
  // exit, the value will be -1).
  std::vector<Int> last_degree_member;

  // Successively increases degree_lower_bound until a degree list member is
  // found (and an index of minimal degree is returned).
  // If 'demand_smallest_index' is true, then the member of the linked list with
  // the smallest index is returned.
  Int FindMinimalIndex(bool demand_smallest_index);

  // Removes the current degree of occurrence with the given index.
  void RemoveDegree(Int index);

  // Adds in an occurrence of the specified index and degree.
  void AddDegree(Int index, Int degree);

  // If the old degree was the same as the current degree, then this routine
  // removes the old occurrence and adds in the new occurrence. Otherwise,
  // it is a no-op.
  void UpdateDegree(Int index, Int degree);
};

// A data structure representing the "quotient graph" interpretation of the
// original graph after eliminating a sequence of vertices. This is the
// primary data structure of the (Approximate) Minimum Degree reordering
// algorithm.
//
// Please see:
//
//   [ADD-96]
//   Patrick R. Amestoy, Timothy A. Davis, and Iain S. Duff,
//   "An Approximate Minimum Degree Ordering Algorithm",
//   SIAM J. Matrix Analysis & Applic., Vol. 17, No. 4, pp. 886--905, 1996.
//
// There appears to be a plethora of publicly available preprints.
//
struct QuotientGraph {
  // The number of vertices in the original graph.
  Int num_original_vertices;

  // The number of vertices that have been eliminated from the original graph.
  Int num_eliminated_vertices;

  // A list of length 'num_original_vertices' of the sizes of each supernode.
  // If index 'i' is not principal, then it is set to zero.
  std::vector<Int> supernode_sizes;

  // A list of length 'num_original_vertices' such that each supernode is
  // traversed from the principal member to the final member in a continuous
  // manner by following the 'next_index' paths.
  //
  // The values are undefined on tail indices of supernodes.
  std::vector<Int> next_index;

  // A list of length 'num_original_vertices' such that, if 'j' is in supernode
  // with principal member 'i', then 'head_index[j]' is 'i'.
  std::vector<Int> head_index;

  // A list of length 'num_original_vertices' such that, if 'i' is principal,
  // then 'tail_index[i]' points to the last index in supernode i (with the
  // ordering defined by the 'next_index' traversal).
  std::vector<Int> tail_index;

  // A list of length 'num_original_vertices' of element nonzero structures.
  // The 'element' index of the list, 'structures[element]', will be created
  // when supernode 'element' is converted from a variable to an element.
  //
  // The structure list also contains non-principal members.
  std::vector<std::vector<Int>> structures;

  // A list of length 'num_original_vertices' of the (unmodified) variable
  // adjacencies of each principal variable. For example, if index 'i' is a
  // principal variable, then 'adjacency_lists[i]' contains the set of neighbor
  // variables for variable i that are not redundant with respect to edges
  // implied by 'structures'.
  //
  // The adjacency lists connect to each individual member of any supernode.
  std::vector<std::vector<Int>> adjacency_lists;

  // A list of length 'num_original_vertices' of the elements adjacent to
  // each principal variable. For example, if index 'i' is a principal
  // variable, then 'element_lists[i]' contains the list of elements adjacent
  // to supervariable 'i'.
  //
  // The element lists only contain the principal member of any supernode.
  std::vector<std::vector<Int>> element_lists;

  // A set of linked lists for keeping track of supervariables of each degree
  // (and, also, a way to provide fast access to a supervariable with
  // minimal degree).
  DegreeLists degree_lists;

  // An optional list of aggressive element absorption pairs: each pair (e, f)
  // consists of the absorbing element, e, and the absorbed element, f.
  std::vector<std::pair<Int, Int>> aggressive_absorptions;

  // An optional list of supervariable merge pairs: each pair (i, j) consists of
  // the absorbing supervariable, i, and the absorbed supervariable, j.
  std::vector<std::pair<Int, Int>> variable_merges;

  // Trivial constructor.
  QuotientGraph();

  // Initializes the quotient graph from a symmetric graph.
  QuotientGraph(const CoordinateGraph& graph);

  // Pretty-prints the QuotientGraph.
  void Print() const;

  // Returns the principal member (if it exists) of the next supernode.
  Int NextSupernode(Int i) const;

  // Returns the principal member (if it exists) of the previous supernode.
  Int PreviousSupernode(Int i) const;

  // Forms the set of members of the supernode with principal variable 'i'.
  std::vector<Int> FormSupernode(Int i) const;

  // Returns a hash of a particular variable.
  std::size_t VariableHash(Int i, VariableHashType hash_type) const;

  // Returns true if supernodes 'i' and 'j' are considered indistinguishable
  // with respect to their quotient graph representation. It is assumed that
  // both supernodes share at least one element as a neighbor.
  //
  // The elimination graph definition (e.g., as given by [ADD-96])
  // involves testing if
  //
  //   Adj_{GElim}(i) \cup {i} = Adj_{GElim}(j) \cup {j}.
  //
  // There is a discussion in [ADD-96] about using
  //
  //   Adj_{GQuotient}(i) \cup {i} = Adj_{GQuotient}(j) \cup {j},
  //
  // but the original George and Liu definition of indistinguishability
  // involved Reach(i) \cup {i} = Reach(j) \cup {j}, where Reach(j) is the
  // union of the adjacencies of node i in the quotient graph *and* its
  // adjacencies that are *through* elements. With this in mind, and the fact
  // that we only query indistinguishability when i and j are known to share
  // an element neighbor, they must be reachable from each other.
  //
  // We therefore test for the equality of the element lists and adjacency
  // lists.
  //
  // TODO(Jack Poulson): Characterize when non-principal members can be in
  // the adjacency lists. In the mean time, we will be conservative and require
  // the non-principal members of the adjacency list to also be equal.
  //
  bool StructuralSupervariablesAreQuotientIndistinguishable(Int i, Int j) const;

  // Initializes 'external_structure_sizes' to a 'num_original_vertices' length
  // vector with all entries equal to -1.
  void InitializeExternalStructureSizes(
    std::vector<Int>* external_structure_sizes) const;

  // An implementation of Algorithm 2 from [ADD-96].
  // On exit, it holds |L_e \ L_p| for all elements e in the element list
  // of a supernode in the structure, L_p.
  //
  // On entry all entries of external_structure_sizes should be less than zero.
  // On exit, all entries of 'external_structure_sizes' corresponding to element
  // indices in the element list of a supernode in the structure L_p should be
  // non-negative and equal to |L_e \ L_p|.
  //
  // If the 'aggressive_absorption' boolean is true, then
  // 'aggressive_absorption_elements' is filled with the elements which should
  // be absorbed.
  void ExternalStructureSizes(
    const std::vector<Int>& supernodal_pivot_structure,
    bool aggressive_absorption,
    std::vector<Int>* external_structure_sizes,
    std::vector<Int>* aggressive_absorption_elements) const;

  // Sets all entries of 'external_structure_sizes' that correspond to an
  // element index in the element list of a supernode in the structure L_p to
  // -1.
  void ResetExternalStructureSizes(
    const std::vector<Int>& supernodal_structure,
    std::vector<Int>* external_structure_sizes) const;

 private:
  // A definition of Ashcraft's hash function (as described in [ADD-96]).
  // Note that only principal members of A_i are incorporated in the hash.
  std::size_t AshcraftVariableHash(Int i) const;

  // An alternative hash that does not explicitly use modular arithmetic and
  // multiplies each index contribution by its position in the adjacency or
  // element list (with the hope of decreasing collisions).
  std::size_t BasicVariableHash(Int i) const;
};

// Pretty-prints an std::vector<T>.
// TODO(Jack Poulson): Find a better location for this utility function.
template<typename T>
void PrintVector(const std::vector<T>& vec, const std::string& msg);

} // namespace quotient

#include "quotient/quotient_graph-impl.hpp"

#endif // ifndef QUOTIENT_QUOTIENT_GRAPH_H_
