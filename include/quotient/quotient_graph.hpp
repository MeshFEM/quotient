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

#include "quotient/buffer.hpp"
#include "quotient/coordinate_graph.hpp"
#include "quotient/degree_and_hash_lists.hpp"
#include "quotient/index_utils.hpp"
#include "quotient/integers.hpp"
#include "quotient/io_utils.hpp"
#include "quotient/macros.hpp"
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
// It is also worth emphasizing that this implementation is similar in many
// respects to that of SuiteSparse's AMD.
//
class QuotientGraph {
 public:
  // Initializes the quotient graph from a symmetric graph.
  QuotientGraph(const CoordinateGraph& graph,
                const MinimumDegreeControl& control);

  // Initializes the quotient graph from the edges of a symmetric graph.
  QuotientGraph(Int num_vertices, const Buffer<GraphEdge>& edges,
                const MinimumDegreeControl& control);

  // Initializes the quotient graph from the entries of a symmetric matrix.
  template <typename Field>
  QuotientGraph(Int num_vertices, const Buffer<MatrixEntry<Field>>& edges,
                const MinimumDegreeControl& control);

  // Returns the number of vertices in the original graph.
  Int NumVertices() const QUOTIENT_NOEXCEPT;

  // Returns the number of vertices that have been eliminated from the graph.
  Int NumEliminatedVertices() const QUOTIENT_NOEXCEPT;

  // Returns the ordered list of eliminated principal variables.
  const std::vector<Int>& EliminationOrder() const QUOTIENT_NOEXCEPT;

  // Fills 'postorder' with the postorder of the assembly tree.
  void ComputePostorder(Buffer<Int>* postorder) const QUOTIENT_NOEXCEPT;

  // Overwrites 'permuted_supernode_sizes' with the sizes of the supernodes
  // in their permuted order.
  void PermutedSupernodeSizes(const Buffer<Int>& inverse_permutation,
                              Buffer<Int>* permuted_supernode_sizes) const
      QUOTIENT_NOEXCEPT;

  // Overwrites 'permuted_member_to_supernode' with a map from the permuted
  // vertex indices to the containing permuted supernode index.
  void PermutedMemberToSupernode(
      const Buffer<Int>& inverse_permutation,
      Buffer<Int>* permuted_member_to_supernode) const QUOTIENT_NOEXCEPT;

  // Overwrites 'permuted_assembly_parents' with the uplinks of the permuted
  // supernodal assembly forest.
  void PermutedAssemblyParents(const Buffer<Int>& permutation,
                               const Buffer<Int>& permuted_member_to_supernode,
                               Buffer<Int>* permuted_assembly_parents) const
      QUOTIENT_NOEXCEPT;

  // Returns the number of times that supervariables have been falsely hashed
  // into the same bucket.
  Int NumHashBucketCollisions() const QUOTIENT_NOEXCEPT;

  // Returns the number of times that supervariables falsely had the same hash
  // value. This is much more rare than hash bucket collision.
  Int NumHashCollisions() const QUOTIENT_NOEXCEPT;

  // Forms the set of members of the supernode with the given principal member.
  Buffer<Int> FormSupernode(Int i) const QUOTIENT_NOEXCEPT;

  // Returns the size of the supernode with the given principal variable.
  Int SupernodeSize(Int i) const QUOTIENT_NOEXCEPT;

  // Returns the element for the given principal member.
  const Buffer<Int> Element(Int i) const QUOTIENT_NOEXCEPT;

  // Returns a copy of the element list of the given principal member.
  Buffer<Int> ElementList(Int i) const QUOTIENT_NOEXCEPT;

  // Finds the next pivot supervariable, forms the corresponding element, and
  // updates the quotient graph. The return value is the principal member of
  // the selected pivot.
  Int FindAndProcessPivot() QUOTIENT_NOEXCEPT;

  // Returns the number of members of the element list of the pivot.
  Int NumPivotElements() const QUOTIENT_NOEXCEPT;

  // Returns the number of degree updates required to process the current pivot.
  Int NumPivotDegreeUpdates() const QUOTIENT_NOEXCEPT;

  // Returns the number of degree updates required to process the current pivot
  // that will involve more than two elements in the element list.
  Int NumPivotDegreeUpdatesWithMultipleElements() const QUOTIENT_NOEXCEPT;

  // Returns the number of nonzeros in the current pivot's columns of the
  // lower-triangular Cholesky factor.
  Int NumPivotCholeskyNonzeros() const QUOTIENT_NOEXCEPT;

  // Returns the number of floating-point operations required for a standard
  // Cholesky factorization to eliminate the current pivot.
  double NumPivotCholeskyFlops() const QUOTIENT_NOEXCEPT;

  // Returns the number of aggressive absorptions that occurred.
  Int NumAggressiveAbsorptions() const QUOTIENT_NOEXCEPT;

  // Returns the number of dense rows that were preprocessed out.
  Int NumDense() const QUOTIENT_NOEXCEPT;

  // Prints the current breakdown of the stage timings. The result will be
  // trivial unless QUOTIENT_ENABLE_TIMERS is defined.
  std::vector<std::pair<std::string, double>> ComponentSeconds() const
      QUOTIENT_NOEXCEPT;

  // This routine should be called after eliminating the non-dense variables,
  // as it updates the assembly forest by combining the dense nodes into a
  // single supernode which becomes the parent of all non-dense roots.
  void CombineDenseNodes() QUOTIENT_NOEXCEPT;

  // Returns an immutable reference to the control structure.
  const MinimumDegreeControl& Control() const QUOTIENT_NOEXCEPT;

 private:
  // Bookkeeping data for the dense supernode (if it exists). It is used at the
  // end of the minimum degree analysis to postprocess the assembly forest so
  // that, if any dense variables were detected, they are all combined and
  // injected as the root of the assembly tree.
  struct DenseSupernode {
    // The number of dense rows that were preprocessed out.
    Int size;

    // If there are any dense nodes, they will eventually be combined into a
    // single supernode with this principal member.
    Int principal_member;
  };

  // An easily-modifiable representation of the supernodes in the quotient
  // graph. Each supernode is maintained as a singly-linked list.
  struct AssemblyForest {
    // A list of length 'num_vertices' of the (signed) sizes of each
    // supernode. If index 'i' is not principal, then it is set to zero; if
    // 'i' is a principal variable, then index 'i' is the size of the supernode:
    // if 'i' is a principal element, the value is negated.
    //
    // Absorbed elements and dense supernode members both are marked via a
    // signed size of '0', but eliminated elements have their assembly parent
    // marked as their parent in the tree, while dense supernode member 'i' has
    // its parent equal to SYMMETRIC_INDEX(i).
    Buffer<Int> signed_supernode_sizes;

    // A (possibly empty) dense supernode.
    DenseSupernode dense_supernode;

    // A list of length 'num_vertices' where index 'e' contains the index of
    // the parent of element 'e' in the elimination forest (if it exists).
    // If element 'e' has no parent, then the value is equal to
    // SYMMETRIC_INDEX(j), where 'j' is the last member of the supernode.
    Buffer<Int> parent_or_tail;
  };

  // A data structure managing an array of length 'num_vertices' which can be
  // used to quickly compute the cardinalities of |L_e \ L_p| for each element
  // e in an element list of a supervariable in the current pivot structure,
  // L_p.
  //
  // The positive values can be quickly 'unset' by increasing a shift such that,
  // in the next iteration, a value is unset if it is less than the shift.
  //
  // It is also used for temporarily flagging variables as within a set.
  struct NodeFlags {
    // A mask of length 'num_vertices' that can be used to quickly compute
    // the cardinalities of |L_e \ L_p| for each element e in an element list of
    // a supervariable in the current pivot structure, L_p.
    //
    // It is also used for temporarily flagging variables as within a set.
    Buffer<Int> flags;

    // The maximum degree that has been constructed so far. Since the external
    // degree updates in each stage will be less than this value, it is used as
    // the amount to increase external_degree_shift_ by at each iteration.
    //
    // TODO(Jack Poulson): This is true for Amestoy and exact degree bounds, but
    // I have not yet checked if it holds for the Gilbert bound.
    Int max_degree;

    // The current datum value for the external degrees (stored within
    // node_flags_). All values should be interpreted relative to the datum
    // value.
    Int shift;

    // The maximum allowable value of the datum until an explicit reset is
    // required.
    Int shift_cap;
  };

  // A packing of the adjacency and element lists, with the element lists
  // occurring first in each member, so that memory allocations are not
  // required during the elimination process. The list is of length
  // 'num_vertices_'. After a supervariable is converted into an element, the
  // metadata (and potentially storage) is repurposed for storing the element
  // structure.
  //
  // Each element list is the set of current children of a principal variable.
  //
  // The adjacency portion of each member contains the (unmodified) variable
  // adjacencies of the principal variable. For example, if index 'i' is a
  // principal variable, then 'adjacency_lists[i]' contains the set of neighbor
  // variables for variable i that are not redundant with respect to edges
  // implied by 'structures'.
  struct QuotientGraphData {
    // The concatentation of the element + adjacency lists of each node.
    Buffer<Int> lists;

    // When index 'i' is a variable, the element list of supervariable 'i' will
    // start at index `element_offsets[i]` of 'lists'. When supervariable 'i'
    // becomes an element, the element will be stored in this location.
    Buffer<Int> element_offsets;

    // The length of the element list of variable i.

    // When index 'i' is a variable, 'element_sizes[i]' will denote the length
    // of the element list. When 'i' is an element, it will denote the number
    // of supervariables in the element.
    Buffer<Int> element_sizes;

    // The length of the variable list of variable i.
    Buffer<Int> adjacency_list_sizes;

    // The position the next element can be stored at.
    Int offset;

    // Returns a mutable pointer to the element list of a given variable.
    Int* ElementList(Int i) QUOTIENT_NOEXCEPT {
      return &lists[element_offsets[i]];
    }

    // Returns an immutable pointer to the element list of a given variable.
    const Int* ElementList(Int i) const QUOTIENT_NOEXCEPT {
      return &lists[element_offsets[i]];
    }

    // Returns a mutable pointer to the structure of a given element.
    Int* ElementData(Int i) QUOTIENT_NOEXCEPT {
      return &lists[element_offsets[i]];
    }

    // Returns an immutable pointer to the structure of a given element.
    const Int* ElementData(Int i) const QUOTIENT_NOEXCEPT {
      return &lists[element_offsets[i]];
    }

    // Returns a mutable pointer to the adjacency list of a given variable.
    Int* AdjacencyList(Int i) QUOTIENT_NOEXCEPT {
      return &lists[element_offsets[i] + element_sizes[i]];
    }

    // Returns an immutable pointer to the adjacency list of a given variable.
    const Int* AdjacencyList(Int i) const QUOTIENT_NOEXCEPT {
      return &lists[element_offsets[i] + element_sizes[i]];
    }

    // Returns an immutable reference to the element list length for
    // supervariable i.
    const Int& ElementListSize(Int i) const QUOTIENT_NOEXCEPT {
      return element_sizes[i];
    }

    // Returns a mutable reference to the element list length for supervariable
    // i.
    Int& ElementListSize(Int i) QUOTIENT_NOEXCEPT { return element_sizes[i]; }

    // Returns an immutable reference to the structure length for an element.
    const Int& ElementSize(Int element) const QUOTIENT_NOEXCEPT {
      return element_sizes[element];
    }

    // Returns a mutable reference to the structure length for an element.
    Int& ElementSize(Int element) QUOTIENT_NOEXCEPT {
      return element_sizes[element];
    }

    // Contiguously pack the still-active elements into 'indices'.
    void PackElements(Int pack_offset, const Int* element_beg,
                      const Int* element_end) QUOTIENT_NOEXCEPT {
      offset = pack_offset;
      for (const Int* iter = element_beg; iter != element_end; ++iter) {
        const Int element = *iter;
        const Int element_size = element_sizes[element];
        Int element_offset = element_offsets[element];
        element_offsets[element] = offset;
        for (Int i = 0; i < element_size; ++i) {
          lists[offset++] = lists[element_offset++];
        }
      }
    }
  };

  // Data structures related to degree lists and hashed supervariables.
  struct DegreesAndHashes {
    // A set of linked lists for keeping track of supervariables of each degree
    // (and, also, a way to provide fast access to a supervariable with
    // minimal degree).
    //
    // Hashes and a singly-linked hash bucket list are stored within the
    // leftovers from removing degree links to principal supervariables in the
    // pivot structure.
    DegreeAndHashLists lists;

    // The number of times that supervariables were falsely placed within the
    // same bucket.
    Int num_bucket_collisions;

    // The number of times that supervariables falsely had the same hash value.
    Int num_collisions;
  };

  // The control structure used to configure the MinimumDegree analysis.
  const MinimumDegreeControl control_;

  // The number of vertices in the original graph.
  Int num_vertices_;

  // The number of vertices that have been eliminated from the original graph.
  Int num_eliminated_vertices_;

  // The principal member of the current pivot.
  Int pivot_;

  // The representation of the current assembly forest.
  AssemblyForest assembly_;

  // The representation of the element lists and adjacencies of the nodes
  // in the quotient graph.
  QuotientGraphData graph_data_;

  // A data structure for quickly maintaining node statuses and degrees.
  NodeFlags node_flags_;

  // Data structures related to degree lists and hashed supervariables.
  DegreesAndHashes degrees_and_hashes_;

  // The number of aggressive absorptions that have occurred.
  Int num_aggressive_absorptions_;

  // The ordered list of principal members of eliminated supernodes.
  std::vector<Int> elimination_order_;

#ifdef QUOTIENT_ENABLE_TIMERS
  // A map from the stage name to the associated timer.
  mutable std::unordered_map<std::string, Timer> timers_;
#endif

  // Initialize the AssemblyForest.
  void InitializeAssemblyForest() QUOTIENT_NOEXCEPT;

  // Converts edge counts for each source into an offset scan and return the
  // number of edges.
  Int ConvertEdgeCountsIntoOffsets() QUOTIENT_NOEXCEPT;

  // Initializes the DegreeAndHashLists data structure.
  void InitializeDegreeAndHashLists() QUOTIENT_NOEXCEPT;

  // Initializes the NodeFlags data structure.
  void InitializeNodeFlags() QUOTIENT_NOEXCEPT;

  // Retrieve a variable with minimal (approximate) degree and set it as the
  // active pivot.
  Int GetNextPivot() QUOTIENT_NOEXCEPT;

  // Stores the element for the pivot:
  //
  //   L_p := (A_p \cup (\cup_{e in E_p} L_e)) \ supernode(p).
  //
  // It is assumed that the mask is of length 'num_vertices' and set to all
  // zeros on input.
  void ComputePivotStructure() QUOTIENT_NOEXCEPT;

  // Compute the degree approximations of the supernodes adjacent to the
  // current pivot. Element absorption is performed during this call.
  void ComputeDegreesAndHashes() QUOTIENT_NOEXCEPT;

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
  bool StructuralVariablesAreQuotientIndistinguishable(Int i, Int j) const
      QUOTIENT_NOEXCEPT;

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
  void MergeVariables() QUOTIENT_NOEXCEPT;

  // Performs the final cleanup for the processing of the pivot element.
  void FinalizePivot() QUOTIENT_NOEXCEPT;

  // An implementation of Algorithm 2 from [ADD-96].
  // On exit, it holds |L_e \ L_p| for all elements e in the element list
  // of a supernode in the structure, L_p.
  //
  // On entry all entries of node_flags_ should be less than the external
  // element size shift.
  //
  // On exit, all entries of 'node_flags_' corresponding to element indices
  // in the element list of a supernode in the structure L_p should be,
  // after removing the shift, non-negative and equal to |L_e \ L_p|.
  void ExternalDegrees() QUOTIENT_NOEXCEPT;

  // Sets all entries of 'node_flags_' that correspond to an element in the
  // element list of a supernode in the pivot structure, L_p.
  void ResetExternalDegrees() QUOTIENT_NOEXCEPT;

  // Uses the parents_ links for the assembly tree to contiguously fill a
  // subtree of the post-order rooted at 'index' using the iterator.
  Int* PreorderTree(Int index, const Buffer<Int>& nonprincipal_members,
                    const Buffer<Int>& nonprincipal_offsets,
                    const Buffer<Int>& children,
                    const Buffer<Int>& child_offsets,
                    Int* iter) const QUOTIENT_NOEXCEPT;

  // A definition of Ashcraft's hash function (as described in [ADD-96]).
  UInt AshcraftVariableHash(Int i) const QUOTIENT_NOEXCEPT;

  // An alternative hash that does not explicitly use modular arithmetic and
  // multiplies each index contribution by its position in the adjacency or
  // element list (with the hope of decreasing collisions).
  UInt BasicVariableHash(Int i) const QUOTIENT_NOEXCEPT;

  // Accumulates the sum of the supernode sizes in the adjacency list and the
  // hash of their indices. While doing so, the non-principal and redundant
  // members are removed (with the remainder packed at the given index).
  void PackCountAndHashAdjacencies(Int i, Int num_elements, Int* degree,
                                   UInt* hash) QUOTIENT_NOEXCEPT;

  // Computes the exact degree of supernode, say, i, using a short-cut of
  // Eq. (2) of [ADD-96] meant for the case where there is only one member of
  // the element list.
  //   d_i = |A_i \ supernode(i)| + |L_p \ supernode(i)|.
  //
  // NOTE: It is assumed that this supervariable is in the pivot structure.
  std::pair<Int, UInt> ExactEmptyDegreeAndHash(Int i) QUOTIENT_NOEXCEPT;

  // Computes the exact degree of supernode i using a short-cut of Eq. (2) of
  // [ADD-96] meant for the case where there are two members of the element
  // list.
  //   d_i = |A_i \ supernode(i)| + |L_p \ supernode(i)| + |L_e \ L_p|.
  //
  // NOTE: It is assumed that this supervariable is in the pivot structure.
  std::pair<Int, UInt> ExactSingleDegreeAndHash(Int i) QUOTIENT_NOEXCEPT;

  // Computes the exact degree of supernode i using Eq. (2) of [ADD-96] in the
  // case of arbitrary members in element_lists[i].
  //   d_i = |A_i \ supernode(i)| + |(\cup_{e in E_i) L_e) \ supernode(i)|.
  //
  // NOTE: It is assumed that this supervariable is in the pivot structure.
  std::pair<Int, UInt> ExactGenericDegreeAndHash(Int i) QUOTIENT_NOEXCEPT;

  // Updates the exact degrees and hashes of the principal members of the
  // supernodes in the pivot structure.
  void ExactDegreesAndHashes() QUOTIENT_NOEXCEPT;

  // Updates the Amestoy degrees and hashes of the principal members of the
  // supernodes in the pivot structure.
  void AmestoyDegreesAndHashes() QUOTIENT_NOEXCEPT;

  // Returns the degree and hash of supernode i in the current pivot structure.
  std::pair<Int, UInt> GilbertDegreeAndHash(Int i) QUOTIENT_NOEXCEPT;

  // Updates the Gilbert degrees and hashes of the principal members of the
  // supernodes in the pivot structure.
  void GilbertDegreesAndHashes() QUOTIENT_NOEXCEPT;

  // Updates the Ashcraft degrees and hashes of the principal members of the
  // supernodes in the pivot structure.
  void AshcraftDegreesAndHashes() QUOTIENT_NOEXCEPT;

  // Inserts the current pivot into the back of the element list of principal
  // variable 'i' by appending the first adjacency to the back of the adjacency
  // list then replacing the first adjacency with the pivot.
  void InsertPivotElement(Int i) QUOTIENT_NOEXCEPT;
};

}  // namespace quotient

#include "quotient/quotient_graph-impl.hpp"

#endif  // ifndef QUOTIENT_QUOTIENT_GRAPH_H_
