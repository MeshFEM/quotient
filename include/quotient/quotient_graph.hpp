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
#include "quotient/degree_lists.hpp"
#include "quotient/minimum_degree_control.hpp"

namespace quotient {

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
class QuotientGraph {
 public:
  // Trivial constructor.
  QuotientGraph();

  // Initializes the quotient graph from a symmetric graph.
  QuotientGraph(
      const CoordinateGraph& graph, const MinimumDegreeControl& control);

  // Pretty-prints the QuotientGraph.
  void Print() const;

  // Returns the number of vertices in the original graph.
  Int NumOriginalVertices() const;

  // Returns the number of vertices that have been eliminated from the graph.
  Int NumEliminatedVertices() const;

  // Returns the ordered list of eliminated principal variables.
  const std::vector<Int>& EliminationOrder() const;

  // Returns the number of times that supervariables have been falsely hashed
  // into the same bucket.
  Int NumHashCollisions() const;

  // Forms the set of members of the supernode with principal variable 'i'.
  std::vector<Int> FormSupernode(Int i) const;

  //  If control_.store_structures was true, then this routine overwrites
  // 'eliminated_structures' with the (sorted) structures of the eliminated
  // supernodes, in the order in which they were eliminated.
  void FormEliminatedStructures(
      std::vector<std::vector<Int>>* eliminated_structures) const;

  // Retrieve a variable with minimal (approximate) external degree and set it
  // as the active pivot.
  Int GetNextPivot();

  // Stores the element for the pivot:
  //
  //   L_p := (A_p \cup (\cup_{e in E_p} L_e)) \ supernode(p).
  //
  // It is assumed that the mask is of length 'num_orig_vertices' and set to all
  // zeros on input.
  //
  // The return value is the number of traversed members of the elements in the
  // element list of the pivot that are no longer variables.
  Int ComputePivotStructure();

  // Returns the number of members of the element list of the pivot.
  Int NumPivotElements() const;

  // Returns the number of degree updates required to process the current pivot.
  Int NumPivotDegreeUpdates() const;

  //Returns the number of degree updates required to process the current pivot
  // that will involve more than two elements in the element list.
  Int NumPivotDegreeUpdatesWithMultipleElements() const;

  // Returns the number of nonzeros in the current pivot's columns of the
  // lower-triangular Cholesky factor.
  Int NumPivotCholeskyNonzeros() const;

  // Returns the number of floating-point operations required for a standard
  // Cholesky factorization to eliminate the current pivot.
  double NumPivotCholeskyFlops() const;

  // Update the adjacency lists after computing the supernodal pivot structure.
  void UpdateAdjacencyListsAfterSelectingPivot();

  // Update the element lists after computing the supernodal pivot structure.
  void UpdateElementListsAfterSelectingPivot();

  // Sets the entry mask[i] to zero for each i in indices.
  void UnflagPivotStructure();

  // Mark the pivot elements in the mask.
  void FlagPivotElementList();

  // Unmark the pivot elements in the mask.
  void UnflagPivotElementList();

  // Perform any aggressive absorptions and clear the members of the element
  // lists of the absorbed elements from the pivot mask.
  void AggressiveAbsorption(
      const std::vector<Int>& aggressive_absorption_elements);

  // Returns (an approximation of) the external degree of a given supervariable.
  Int ExternalDegree(Int i) const;

  // Compute the external degree approximations of the supernodes adjacent to
  // the current pivot.
  void ComputeExternalDegrees(std::vector<Int>* external_degrees);

  // Insert the new external degrees in to the degree lists.
  void UpdateExternalDegrees(const std::vector<Int>& external_degrees);

  // Returns a hash of a particular variable.
  std::size_t VariableHash(Int i, VariableHashType hash_type) const;

  // Compute hashes of the supervariables in the pivot structure.
  void ComputeVariableHashes(std::vector<std::size_t>* bucket_keys);

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
  bool StructuralSupervariablesAreQuotientIndistinguishable(Int i, Int j) const;

  // Detects and merges pairs of supervariables in the pivot structure who are
  // indistinguishable with respect to the quotient graph.
  //
  // While the supernodal merges will potentially shrink the supernodal
  // adjacency lists (and thus change the associated Ashcraft hash of
  // variables), if two variables are indistinguishable, their cached bucket
  // might be wrong, but they would be wrong together.
  //
  // The test for indistinguishability does not depend upon the variable
  // supernodal structure and is thus invariant to supervariable merges.
  void MergeVariables(const std::vector<std::size_t>& bucket_keys);

  // Converts the 'pivot' (super)variable into an element.
  void ConvertPivotIntoElement();

  // An implementation of Algorithm 2 from [ADD-96].
  // On exit, it holds |L_e \ L_p| for all elements e in the element list
  // of a supernode in the structure, L_p.
  //
  // On entry all entries of external_element_sizes should be less than zero.
  // On exit, all entries of 'external_element_sizes' corresponding to element
  // indices in the element list of a supernode in the structure L_p should be
  // non-negative and equal to |L_e \ L_p|.
  //
  // If the 'aggressive_absorption' boolean is true, then
  // 'aggressive_absorption_elements' is filled with the elements which should
  // be absorbed.
  void RecomputeExternalElementSizes(
    std::vector<Int>* aggressive_absorption_elements);

  // Sets all entries of 'external_element_sizes' that correspond to an
  // element index in the element list of a supernode in the structure L_p to
  // -1.
  void ResetExternalElementSizes();

  // Whether or not external element sizes are used for external degree
  // computations.
  bool UsingExternalElementSizes() const;

  // Returns the list of supervariable merge pairs
  // (if control.store_variable_merges was true): each pair (i, j) consists of
  // the absorbing supervariable, i, and the absorbed supervariable, j.
  const std::vector<std::pair<Int, Int>>& VariableMerges() const;

  // Returns an optional list of aggressive element absorption pairs: each
  // pair (e, f) consists of the absorbing element, e, and the absorbed
  // element, f. This will only return a non-empty list if
  // control.store_aggressive_absorptions was true.
  const std::vector<std::pair<Int, Int>>& AggressiveAbsorptions() const;

 private:
  // A data structure for representing the relevant information of a merge of
  // one supervariable (the 'absorbed_index') into another
  // (the 'primary_index').
  struct VariableMergeInfo {
    // The principal index of the supervariable that is increasing in size.
    Int primary_index;

    // The principal index of the supervariable being absorbed.
    Int absorbed_index;

    // The number of members of the absorbed supervarible.
    Int absorbed_size;

    VariableMergeInfo(
        Int primary_index_value,
        Int absorbed_index_value,
        Int absorbed_size_value) :
    primary_index(primary_index_value),
    absorbed_index(absorbed_index_value),
    absorbed_size(absorbed_size_value) {}
  };

  // The number of vertices in the original graph.
  Int num_original_vertices_;

  // The number of vertices that have been eliminated from the original graph.
  Int num_eliminated_vertices_;

  // The control structure used to configure the MinimumDegree analysis.
  MinimumDegreeControl control_; 

  // Whether or not external element sizes will be maintained throughout the
  // minimum degree analysis.
  bool using_external_element_sizes_;

  // The principal member of the current pivot.
  Int pivot_;

  // A list of length 'num_original_vertices' of the (signed) sizes of each
  // supernode. If index 'i' is not principal, then it is set to zero; if
  // 'i' is a principal variable, then 'supernode_sizes[i]' is the size of the
  // supernode; if 'i' is a principal element, the value is negated.
  std::vector<Int> supernode_sizes_;

  // The ordered list of principal members of eliminated supernodes.
  std::vector<Int> elimination_order_;

  // A list of length 'num_original_vertices' such that each supernode is
  // traversed from the principal member to the final member in a continuous
  // manner by following the 'next_index' paths.
  //
  // The values are undefined on tail indices of supernodes.
  std::vector<Int> next_index_;

  // A list of length 'num_original_vertices' such that, if 'i' is principal,
  // then 'tail_index[i]' points to the last index in supernode i (with the
  // ordering defined by the 'next_index' traversal).
  std::vector<Int> tail_index_;

  // A list of length 'num_original_vertices' of the (unmodified) variable
  // adjacencies of each principal variable. For example, if index 'i' is a
  // principal variable, then 'adjacency_lists[i]' contains the set of neighbor
  // variables for variable i that are not redundant with respect to edges
  // implied by 'structures'.
  //
  // The adjacency lists connect to each individual member of any supernode.
  std::vector<std::vector<Int>> adjacency_lists_;

  // A list of length 'num_original_vertices' of the elements adjacent to
  // each principal variable. For example, if index 'i' is a principal
  // variable, then 'element_lists[i]' contains the list of elements adjacent
  // to supervariable 'i'.
  //
  // The element lists only contain the principal member of any supernode.
  std::vector<std::vector<Int>> element_lists_;

  // A copy of the pivot element list that can be used to update 'element_sizes'
  // when the pivot is converted to an element. This is only needed when
  // aggressive absorption is activated, as the element list will no longer
  // correspond to the set that needs to be decremented by the pivot supernode
  // size.
  std::vector<Int> original_pivot_element_list_;

  // A set of linked lists for keeping track of supervariables of each degree
  // (and, also, a way to provide fast access to a supervariable with
  // minimal degree).
  DegreeLists degree_lists_;

  // An (optional) list of length 'num_original_vertices' of element nonzero
  // structures. The 'element' index of the list, 'structures[element]', will
  // be created when supernode 'element' is converted from a variable to an
  // element.
  //
  // The structure list also contains non-principal members.
  std::vector<std::vector<Int>> structures_;

  // A list of length 'num_original_vertices' of elements (lists of principal
  // variables in the nonzero pattern). The 'e' index of the list,
  // 'elements[e]', will be created when supernode 'e' is converted from a
  // variable to an element.
  //
  // The structure list also contains non-principal members.
  std::vector<std::vector<Int>> elements_;

  // The total number of variables (including nonprincipal) in each element.
  std::vector<Int> element_sizes_;

  // A mask of length 'num_orig_vertices' that can be used to quickly compute
  // the cardinalities of |L_e \ L_p| for each element e in an element list of
  // a supervariable in the current pivot structure, L_p.
  std::vector<Int> external_element_sizes_;

  // A mask of length 'num_original_vertices' that is 1 in index 'i' if and only
  // if 'i' is a member of the current pivot's structure. All other entries
  // will be zero.
  std::vector<int> pivot_mask_;

  // A mask of length 'num_original_vertices' that used within exact external
  // degree computations to perform set unions. It is only created if exact
  // degree computations were requested, and it must be set to all zeros before
  // and after each call to ExternalDegree.
  mutable std::vector<int> exact_degree_mask_;

  // A set of buckets for each hash value (modulo num_original_vertices) of the
  // supervariables.
  std::vector<std::vector<Int>> buckets_;

  // The number of times that supervariables were falsely placed within the
  // same bucket.
  Int num_hash_collisions_;

  // An optional list of supervariable merge pairs: each pair (i, j) consists of
  // the absorbing supervariable, i, and the absorbed supervariable, j.
  std::vector<std::pair<Int, Int>> variable_merges_;

  // An optional list of aggressive element absorption pairs: each pair (e, f)
  // consists of the absorbing element, e, and the absorbed element, f.
  std::vector<std::pair<Int, Int>> aggressive_absorptions_;

  // A definition of Ashcraft's hash function (as described in [ADD-96]).
  // Note that only principal members of A_i are incorporated in the hash.
  std::size_t AshcraftVariableHash(Int i) const;

  // An alternative hash that does not explicitly use modular arithmetic and
  // multiplies each index contribution by its position in the adjacency or
  // element list (with the hope of decreasing collisions).
  std::size_t BasicVariableHash(Int i) const;

  // Computes the exact external degree of supernode i using a short-cut of
  // Eq. (2) of [ADD-96] meant for the case where there are no members in the
  // element list.
  //   d_i = |A_i \ supernode(i)|.
  Int ExactEmptyExternalDegree(Int i) const;

  // Computes the exact external degree of supernode i using a short-cut of
  // Eq. (2) of [ADD-96] meant for the case where there is only one member of
  // the element list.
  //   d_i = |A_i \ supernode(i)| + |L_p \ supernode(i)|.
  //
  // NOTE: It is assumed that 'i' is a principal member of L_p.
  Int ExactSingleExternalDegree(Int i) const;

  // Computes the exact external degree of supernode i using a short-cut of
  // Eq. (2) of [ADD-96] meant for the case where there are two members of the
  // element list.
  //   d_i = |A_i \ supernode(i)| + |L_p \ supernode(i)| + |L_e \ L_p|.
  //
  // NOTE: It is assumed that 'i' is a principal member of L_p.
  Int ExactDoubleExternalDegree(Int i) const;

  // Computes the exact external degree of supernode i using Eq. (2) of
  // [ADD-96] in the case of arbitrary members in element_lists[i].
  //   d_i = |A_i \ supernode(i)| + |(\cup_{e in E_i) L_e) \ supernode(i)|.
  Int ExactGenericExternalDegree(Int i) const;

  // Computes the exact external degree of supernode i using Eq. (2) of
  // [ADD-96].
  //   d_i = |A_i \ supernode(i)| + |(\cup_{e in E_i) L_e) \ supernode(i)|.
  Int ExactExternalDegree(Int i) const;

  // Computes an approximation of the external degree of supernode i using
  // Eq. (4) of [ADD-96].
  //
  // This routine is made slightly more accurate by subtracting the size of
  // supernode(i) from bound0.
  //
  // NOTE: It is assumed that 'i' is a principal member of L_p.
  Int AmestoyExternalDegree(Int i) const;

  // Returns the external degree approximation of Gilbert, Moler, and Schreiber,
  //   \hat{d_i} = |A_i \ supernode(i)|  + \sum_{e in E_i} |L_e \ supernode(i)|.
  //
  // We slightly modify this formula to ensure that the estimated degree is
  // at most
  //   (num_original_vertices - num_eliminated_vertices) - size(supernode(i)).
  Int GilbertExternalDegree(Int i) const;

  // Returns the external degree approximation of Ashcraft, Eisenstat, and
  // Lucas:
  //   \tilde{d_i} = d_i if |E_i| = 2, \hat{d_i} otherwise.
  Int AshcraftExternalDegree(Int i) const;
};

// Pretty-prints an std::vector<T>.
// TODO(Jack Poulson): Find a better location for this utility function.
template<typename T>
void PrintVector(const std::vector<T>& vec, const std::string& msg);

} // namespace quotient

#include "quotient/quotient_graph-impl.hpp"

#endif // ifndef QUOTIENT_QUOTIENT_GRAPH_H_
