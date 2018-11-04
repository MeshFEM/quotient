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
#include "quotient/hash_lists.hpp"
#include "quotient/minimum_degree_control.hpp"
#include "quotient/timer.hpp"

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

  // Fills 'postorder' with the postorder of the assembly tree.
  void ComputePostorder(std::vector<Int>* postorder) const;

  // Returns a reference to the list containing the parent of each supernode.
  const std::vector<Int>& Parents() const;

  // Returns the number of times that supervariables have been falsely hashed
  // into the same bucket.
  Int NumHashBucketCollisions() const;

  // Returns the number of times that supervariables falsely had the same hash
  // value. This is much more rare than hash bucket collision.
  Int NumHashCollisions() const;

  // Forms the set of members of the supernode with the given principal member.
  std::vector<Int> FormSupernode(Int i) const;

  // Returns a reference to the element for the given principal member.
  const std::vector<Int>& Element(Int i) const;

  // Returns a reference to the element list of the given principal member.
  std::vector<Int> ElementList(Int i) const;

  //  If control_.store_structures was true, then this routine overwrites
  // 'eliminated_structures' with the (sorted) structures of the eliminated
  // supernodes, in the order in which they were eliminated.
  void FormEliminatedStructures(
      std::vector<std::vector<Int>>* eliminated_structures) const;

  // Finds the next pivot supervariable, forms the corresponding element, and
  // updates the quotient graph. The return value is the principal member of
  // the selected pivot.
  Int FindAndProcessPivot();

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

  // Returns the number of aggressive absorptions that occurred.
  Int NumAggressiveAbsorptions() const;

  // Prints the current breakdown of the stage timings. The result will be
  // trivial unless QUOTIENT_ENABLE_TIMERS is defined.
  std::vector<std::pair<std::string, double>> ComponentSeconds() const;

 private:
  // The number of vertices in the original graph.
  Int num_original_vertices_;

  // The number of vertices that have been eliminated from the original graph.
  Int num_eliminated_vertices_;

  // The control structure used to configure the MinimumDegree analysis.
  const MinimumDegreeControl control_;

  // The principal member of the current pivot.
  Int pivot_;

  // A list of length 'num_original_vertices' of the (signed) sizes of each
  // supernode. If index 'i' is not principal, then it is set to zero; if
  // 'i' is a principal variable, then index 'i' is the size of the supernode:
  // if 'i' is a principal element, the value is negated.
  std::vector<Int> signed_supernode_sizes_;

  // The ordered list of principal members of eliminated supernodes.
  std::vector<Int> elimination_order_;

  // A list of length 'num_original_vertices' where index 'e' contains the
  // index of the parent of element 'e' in the elimination forest (if it
  // exists). If element 'e' has no parent, then the value is equal to -1.
  std::vector<Int> parents_;

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

  // A packing of the adjacency and element lists, with the element lists
  // occurring first in each member, so that memory allocations are not
  // required during the elimination process. The list is of length
  // 'num_original_vertices_'.
  //
  // Each element list is the set of current children of a principal variable.
  //
  // The adjacency portion of each member contains the (unmodified) variable
  // adjacencies of the principal variable. For example, if index 'i' is a
  // principal variable, then 'adjacency_lists[i]' contains the set of neighbor
  // variables for variable i that are not redundant with respect to edges
  // implied by 'structures'.
  std::vector<Int> element_and_adjacency_lists_;
  std::vector<Int> element_list_offsets_; 
  std::vector<Int> element_list_sizes_;
  std::vector<Int> adjacency_list_sizes_;

  // A set of linked lists for keeping track of supervariables of each degree
  // (and, also, a way to provide fast access to a supervariable with
  // minimal degree).
  DegreeLists degree_lists_;

  // The maximum degree that has been constructed so far. Since the external
  // degree updates in each stage will be less than this value, it is used as
  // the amount to increase external_degree_shift_ by at each iteration.
  Int max_degree_;

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

  // The current datum value for the external_degrees_. All values should
  // be interpreted relative to the datum value.
  Int external_degree_shift_;

  // The maximum allowable value of the datum until an explicit reset is
  // required.
  Int max_shift_value_;

  // A mask of length 'num_orig_vertices' that can be used to quickly compute
  // the cardinalities of |L_e \ L_p| for each element e in an element list of
  // a supervariable in the current pivot structure, L_p.
  //
  // It is also used for temporarily flagging variables as within a set.
  std::vector<Int> node_flags_;

  // An array of single-linked lists for hash buckets for the supervariables.
  HashLists hash_lists_;

  // The number of times that supervariables were falsely placed within the
  // same bucket.
  Int num_hash_bucket_collisions_;

  // The number of times that supervariables falsely had the same hash value.
  Int num_hash_collisions_;

  // The number of aggressive absorptions that have occurred.
  Int num_aggressive_absorptions_;

#ifdef QUOTIENT_ENABLE_TIMERS
  // A map from the stage name to the associated timer.
  mutable std::unordered_map<std::string, Timer> timers_;
#endif

  // Retrieve a variable with minimal (approximate) external degree and set it
  // as the active pivot.
  Int GetNextPivot();

  // Stores the element for the pivot:
  //
  //   L_p := (A_p \cup (\cup_{e in E_p} L_e)) \ supernode(p).
  //
  // It is assumed that the mask is of length 'num_orig_vertices' and set to all
  // zeros on input.
  void ComputePivotStructure();

  // Compute the degree approximations of the supernodes adjacent to the
  // current pivot. Element absorption is performed during this call.
  void ComputeDegreesAndHashes();

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
  bool StructuralVariablesAreQuotientIndistinguishable(Int i, Int j) const;

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
  void MergeVariables();

  // Performs the final cleanup for the processing of the pivot element.
  void FinalizePivot();

  // An implementation of Algorithm 2 from [ADD-96].
  // On exit, it holds |L_e \ L_p| for all elements e in the element list
  // of a supernode in the structure, L_p.
  //
  // On entry all entries of external_degrees_ should be less than the
  // external element size shift.
  //
  // On exit, all entries of 'node_flags_' corresponding to element indices
  // in the element list of a supernode in the structure L_p should be,
  // after removing the shift, non-negative and equal to |L_e \ L_p|.
  void ExternalDegrees();

  // Sets all entries of 'node_flags_' that correspond to an element in the
  // element list of a supernode in the pivot structure, L_p.
  void ResetExternalDegrees();

  // Appends the supernode with the given principal member and length into
  // a given vector.
  void AppendSupernode(Int i, Int supernode_size, std::vector<Int>* vec) const;

  // Uses the parents_ links for the assembly tree to contiguously fill a
  // subtree of the post-order rooted at 'index' using the iterator.
  std::vector<Int>::iterator PreorderTree(
      Int index,
      const std::vector<Int>& children,
      const std::vector<Int>& child_offsets,
      std::vector<Int>::iterator iter) const;

  // A definition of Ashcraft's hash function (as described in [ADD-96]).
  UInt AshcraftVariableHash(Int i) const;

  // An alternative hash that does not explicitly use modular arithmetic and
  // multiplies each index contribution by its position in the adjacency or
  // element list (with the hope of decreasing collisions).
  UInt BasicVariableHash(Int i) const;

  // Accumulates the sum of the supernode sizes in the adjacency list and the
  // hash of their indices. While doing so, the non-principal and redundant
  // members are removed (with the remainder packed at the given index).
  void PackCountAndHashAdjacencies(
      Int i, Int num_elements, Int* degree, UInt* hash);

  // Computes the exact external degree of supernode, say, i, using a short-cut
  // of Eq. (2) of [ADD-96] meant for the case where there is only one member of
  // the element list.
  //   d_i = |A_i \ supernode(i)| + |L_p \ supernode(i)|.
  //
  // NOTE: It is assumed that this supervariable is in the pivot structure.
  std::pair<Int, UInt> ExactEmptyDegreeAndHash(Int i);

  // Computes the exact external degree of supernode i using a short-cut of
  // Eq. (2) of [ADD-96] meant for the case where there are two members of the
  // element list.
  //   d_i = |A_i \ supernode(i)| + |L_p \ supernode(i)| + |L_e \ L_p|.
  //
  // NOTE: It is assumed that this supervariable is in the pivot structure.
  std::pair<Int, UInt> ExactSingleDegreeAndHash(Int i);

  // Computes the exact external degree of supernode i using Eq. (2) of
  // [ADD-96] in the case of arbitrary members in element_lists[i].
  //   d_i = |A_i \ supernode(i)| + |(\cup_{e in E_i) L_e) \ supernode(i)|.
  //
  // NOTE: It is assumed that this supervariable is in the pivot structure.
  std::pair<Int, UInt> ExactGenericDegreeAndHash(Int i);

  // Updates the exact degrees and hashes of the principal members of the
  // supernodes in the pivot structure.
  void ExactDegreesAndHashes();

  // Updates the Amestoy degrees and hashes of the principal members of the
  // supernodes in the pivot structure.
  void AmestoyDegreesAndHashes();

  // Returns the degree and hash of supernode i in the current pivot structure.
  std::pair<Int, UInt> GilbertDegreeAndHash(Int i);

  // Updates the Gilbert degrees and hashes of the principal members of the
  // supernodes in the pivot structure.
  void GilbertDegreesAndHashes();

  // Updates the Ashcraft degrees and hashes of the principal members of the
  // supernodes in the pivot structure.
  void AshcraftDegreesAndHashes();

  // Inserts the current pivot into the back of the element list of principal
  // variable 'i' by appending the first adjacency to the back of the adjacency
  // list then replacing the first adjacency with the pivot.
  void InsertPivotElement(Int i);
};

// Pretty-prints an std::vector<T>.
// TODO(Jack Poulson): Find a better location for this utility function.
template<typename T>
void PrintVector(const std::vector<T>& vec, const std::string& msg);

} // namespace quotient

#include "quotient/quotient_graph-impl.hpp"

#endif // ifndef QUOTIENT_QUOTIENT_GRAPH_H_
