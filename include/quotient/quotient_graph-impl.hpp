/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_QUOTIENT_GRAPH_IMPL_H_
#define QUOTIENT_QUOTIENT_GRAPH_IMPL_H_

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>

#include "quotient/coordinate_graph.hpp"
#include "quotient/integers.hpp"
#include "quotient/macros.hpp"

#include "quotient/quotient_graph.hpp"

namespace quotient {

#ifdef QUOTIENT_ENABLE_TIMERS
#define QUOTIENT_START_TIMER(timer, name) timer[name].Start()
#define QUOTIENT_STOP_TIMER(timer, name) timer[name].Stop()
#else
#define QUOTIENT_START_TIMER(timer, name)
#define QUOTIENT_STOP_TIMER(timer, name)
#endif

#ifdef QUOTIENT_ENABLE_TIMERS
static constexpr char kSetup[] = "Setup";
static constexpr char kComputePivotStructure[] = "ComputePivotStructure";
static constexpr char kExternalDegrees[] = "ExternalDegrees";
static constexpr char kComputeDegrees[] = "ComputeDegrees";
static constexpr char kMergeVariables[] = "MergeVariables";
static constexpr char kFinalizePivot[] = "FinalizePivot";
#endif

inline QuotientGraph::QuotientGraph(const CoordinateGraph& graph,
                                    const MinimumDegreeControl& control)
    : control_(control),
      num_vertices_(graph.NumSources()),
      num_eliminated_vertices_(0),
      num_aggressive_absorptions_(0) {
  QUOTIENT_START_TIMER(timers_, kSetup);

  // Initialize the assembly tree.
  assembly_.signed_supernode_sizes.Resize(num_vertices_, 1);
  assembly_.dense_supernode.size = 0;
  assembly_.dense_supernode.principal_member = -1;
  assembly_.parent_or_tail.Resize(num_vertices_);
  for (Int i = 0; i < num_vertices_; ++i) {
    assembly_.parent_or_tail[i] = SYMMETRIC_INDEX(i);
  }

  elimination_order_.reserve(num_vertices_);

  // Initialize the counts for the adjacency lists.
  edges_.element_list_offsets.Resize(num_vertices_, 0);
  edges_.element_list_sizes.Resize(num_vertices_, 0);
  edges_.adjacency_list_sizes.Resize(num_vertices_, 0);
  const Buffer<GraphEdge>& edges = graph.Edges();
  for (Int source = 0; source < num_vertices_; ++source) {
    const Int source_edge_offset = graph.SourceEdgeOffset(source);
    const Int next_source_edge_offset = graph.SourceEdgeOffset(source + 1);
    for (Int k = source_edge_offset; k < next_source_edge_offset; ++k) {
      const GraphEdge& edge = edges[k];
      if (edge.second != source) {
        ++edges_.adjacency_list_sizes[source];
      }
    }
  }

  // Every row with at least this many non-diagonal nonzeros will be treated
  // as dense and moved to the back of the postordering.
  const Int dense_threshold =
      std::max(1.f * control_.min_dense_threshold,
               control_.dense_sqrt_multiple * std::sqrt(1.f * num_vertices_));

  // Convert the counts into an offset scan.
  Int num_edges = 0;
  for (Int source = 0; source < num_vertices_; ++source) {
    edges_.element_list_offsets[source] = num_edges;
    if (edges_.adjacency_list_sizes[source] >= dense_threshold) {
      if (!assembly_.dense_supernode.size) {
        assembly_.dense_supernode.principal_member = source;
      }
      ++assembly_.dense_supernode.size;
      ++num_eliminated_vertices_;
      // We will denote this source as dense by setting its supernode size to
      // 0 and keeping its parent as -1. A merged variable also has a supernode
      // size of 0, but its parent is a valid vertex index.
      assembly_.signed_supernode_sizes[source] = 0;
      edges_.adjacency_list_sizes[source] = 0;
    }
    num_edges += edges_.adjacency_list_sizes[source];
  }
#ifdef QUOTIENT_DEBUG
  if (assembly_.dense_supernode.size) {
    std::cout << "Eliminated " << assembly_.dense_supernode.size
              << " dense rows." << std::endl;
  }
#endif

  // Pack the edges.
  edges_.lists.Resize(num_edges);
  num_edges = 0;
  for (Int source = 0; source < num_vertices_; ++source) {
    if (!assembly_.signed_supernode_sizes[source]) {
      // Skip the dense row.
      continue;
    }
    const Int source_edge_offset = graph.SourceEdgeOffset(source);
    const Int next_source_edge_offset = graph.SourceEdgeOffset(source + 1);
    for (Int k = source_edge_offset; k < next_source_edge_offset; ++k) {
      const GraphEdge& edge = edges[k];
      if (edge.second != source) {
        edges_.lists[num_edges++] = edge.second;
      }
    }
  }

  // Initialize the degree lists.
  degree_lists_.degrees.Resize(num_vertices_, 0);
  degree_lists_.heads.Resize(num_vertices_ - 1, -1);
  degree_lists_.next_member.Resize(num_vertices_, -1);
  degree_lists_.last_member.Resize(num_vertices_, -1);
  for (Int source = 0; source < num_vertices_; ++source) {
    if (!assembly_.signed_supernode_sizes[source]) {
      // Skip the dense row.
      continue;
    }
    const Int degree = edges_.adjacency_list_sizes[source];
    degree_lists_.AddDegree(source, degree);
  }

  // Initialize the hash lists.
  hash_info_.num_collisions = 0;
  hash_info_.num_bucket_collisions = 0;
  hash_info_.lists.buckets.Resize(num_vertices_);
  hash_info_.lists.hashes.Resize(num_vertices_);
  hash_info_.lists.heads.Resize(num_vertices_, -1);
  hash_info_.lists.next_member.Resize(num_vertices_, -1);

  // Trivially initialize the lower-triangular nonzero structures.
  elements_.Resize(num_vertices_);

  // The absorbed elements will be maintained as 0, the unabsorbed will be
  // initialized as 1, and the shift will be initialized as 2.
  node_flags_.shift = 2;
  node_flags_.flags.Resize(num_vertices_, 1);
  node_flags_.max_degree = 0;
  node_flags_.shift_cap = std::numeric_limits<Int>::max() - num_vertices_;

  QUOTIENT_STOP_TIMER(timers_, kSetup);
}

inline const std::vector<Int>& QuotientGraph::EliminationOrder() const
    QUOTIENT_NOEXCEPT {
  return elimination_order_;
}

inline Int QuotientGraph::FindAndProcessPivot() QUOTIENT_NOEXCEPT {
  // Get the next pivot.
  GetNextPivot();

  // Compute the structure of this pivot block.
  ComputePivotStructure();

  // Compute the external structure cardinalities, |L_e \ L_p|, of all
  // elements e in an element list of a supernode in L_p.
  ExternalDegrees();

  // Store the external degrees of all supervariables in the pivot structure.
  ComputeDegreesAndHashes();

  if (control_.allow_supernodes) {
    // Merge any equivalent supernodes by explicitly checking for equality
    // between pairs that are in the same hash bucket.
    MergeVariables();
  }

  FinalizePivot();

  return pivot_;
}

inline Int QuotientGraph::GetNextPivot() QUOTIENT_NOEXCEPT {
  pivot_ = degree_lists_.FindMinimalIndex(control_.force_minimal_pivot_indices);
  return pivot_;
}

inline Int QuotientGraph::NumAggressiveAbsorptions() const QUOTIENT_NOEXCEPT {
  return num_aggressive_absorptions_;
}

inline Int QuotientGraph::NumDense() const QUOTIENT_NOEXCEPT {
  return assembly_.dense_supernode.size;
}

inline Int QuotientGraph::NumVertices() const QUOTIENT_NOEXCEPT {
  return num_vertices_;
}

inline Int QuotientGraph::NumEliminatedVertices() const QUOTIENT_NOEXCEPT {
  return num_eliminated_vertices_;
}

inline Int QuotientGraph::NumHashBucketCollisions() const QUOTIENT_NOEXCEPT {
  return hash_info_.num_bucket_collisions;
}

inline Int QuotientGraph::NumHashCollisions() const QUOTIENT_NOEXCEPT {
  return hash_info_.num_collisions;
}

inline Int QuotientGraph::SupernodeSize(Int i) const QUOTIENT_NOEXCEPT {
  if (assembly_.signed_supernode_sizes[i]) {
    // This is a traditional supernode (but it might be eliminated).
    return std::abs(assembly_.signed_supernode_sizes[i]);
  } else if (assembly_.parent_or_tail[i] < 0) {
    // This is a dense node.
    return 1;
  } else {
    // This is a merged node.
    return 0;
  }
}

inline const Buffer<Int>& QuotientGraph::Element(Int i) const
    QUOTIENT_NOEXCEPT {
  return elements_[i];
}

inline Buffer<Int> QuotientGraph::ElementList(Int i) const QUOTIENT_NOEXCEPT {
  const Int num_elements = edges_.element_list_sizes[i];
  const Int element_list_beg = edges_.element_list_offsets[i];
  const Int element_list_end = element_list_beg + num_elements;

  Buffer<Int> element_list(num_elements);
  Int counter = 0;
  for (Int k = element_list_beg; k < element_list_end; ++k) {
    element_list[counter++] = edges_.lists[k];
  }

  return element_list;
}

inline void QuotientGraph::ComputePivotStructure() QUOTIENT_NOEXCEPT {
  QUOTIENT_START_TIMER(timers_, kComputePivotStructure);
  QUOTIENT_ASSERT(elements_[pivot_].Empty(), "Chose a pivot more than once.");
  const Int element_list_beg = edges_.element_list_offsets[pivot_];
  const Int element_list_end =
      element_list_beg + edges_.element_list_sizes[pivot_];
  const Int adjacency_list_end =
      element_list_end + edges_.adjacency_list_sizes[pivot_];
  const Int pivot_supernode_size = assembly_.signed_supernode_sizes[pivot_];
  Buffer<Int>& pivot_element = elements_[pivot_];

  // Allocate space for the element using an upper-bound on the size
  // (note that, because of supernodes, this is *not* the degree).
  Int element_size_bound = edges_.adjacency_list_sizes[pivot_];
  for (Int k = element_list_beg; k < element_list_end; ++k) {
    const Int element = edges_.lists[k];
    element_size_bound += elements_[element].Size() - 1;
  }
  pivot_element.Resize(element_size_bound);

  // Negate the signed supernode size of the pivot.
  assembly_.signed_supernode_sizes[pivot_] = -pivot_supernode_size;

  // Push the supervariables in the pivot adjacency list into the structure.
  Int pivot_size = 0;
  Int pivot_degree = 0;
  for (Int k = element_list_end; k < adjacency_list_end; ++k) {
    const Int i = edges_.lists[k];
    const Int supernode_size = assembly_.signed_supernode_sizes[i];
    QUOTIENT_ASSERT(supernode_size >= 0,
                    "An element was in the adjacency list.");
    if (!supernode_size) {
      continue;
    }
    QUOTIENT_ASSERT(assembly_.parent_or_tail[i] < 0,
                    "An absorbed element was in the adjacency list.");

    pivot_degree += supernode_size;
    assembly_.signed_supernode_sizes[i] = -supernode_size;
    pivot_element[pivot_size++] = i;
  }
  edges_.adjacency_list_sizes[pivot_] = 0;

  // Push the unique supervariables in the patterns of the pivot element list
  // into the current structure.
  for (Int k = element_list_beg; k < element_list_end; ++k) {
    const Int element = edges_.lists[k];
    QUOTIENT_ASSERT(assembly_.parent_or_tail[element] < 0,
                    "Used an absorbed element in pivot structure.");

    for (const Int& index : elements_[element]) {
      const Int supernode_size = assembly_.signed_supernode_sizes[index];
      // While no eliminated supernodes should appear in unabsorbed elements,
      // we have (temporarily) flipped the signs of the supernode sizes of
      // the members we have already added to this pivot's structure.
      if (supernode_size <= 0) {
        continue;
      }
      QUOTIENT_ASSERT(assembly_.parent_or_tail[index] < 0,
                      "An absorbed element was in an element.");

      pivot_degree += supernode_size;
      assembly_.signed_supernode_sizes[index] = -supernode_size;
      pivot_element[pivot_size++] = index;
    }
#ifdef QUOTIENT_DEBUG
    Int degree = 0;
    for (const Int& index : elements_[element]) {
      const Int supernode_size = -assembly_.signed_supernode_sizes[index];
      QUOTIENT_ASSERT(supernode_size >= 0,
                      "Flipped supernode size was expected to be positive.");
      degree += supernode_size;
    }
    QUOTIENT_ASSERT(degree == degree_lists_.degrees[element],
                    "Degree did not match its cached value.");
#endif

    // Absorb this element into the pivot.
    assembly_.parent_or_tail[element] = pivot_;
    node_flags_.flags[element] = 0;
    elements_[element].Clear();
    edges_.element_list_sizes[element] = 0;
  }

  // Shrink the size of the pivot element to fit its length.
  pivot_element.Resize(pivot_size);

  degree_lists_.degrees[pivot_] = pivot_degree;
  node_flags_.max_degree = std::max(node_flags_.max_degree, pivot_degree);

  QUOTIENT_STOP_TIMER(timers_, kComputePivotStructure);
}

inline Int QuotientGraph::NumPivotElements() const QUOTIENT_NOEXCEPT {
  return edges_.element_list_sizes[pivot_];
}

inline Int QuotientGraph::NumPivotDegreeUpdates() const QUOTIENT_NOEXCEPT {
  return elements_[pivot_].Size();
}

inline Int QuotientGraph::NumPivotDegreeUpdatesWithMultipleElements() const
    QUOTIENT_NOEXCEPT {
  Int num_multi_updates = 0;
  for (const Int& i : elements_[pivot_]) {
    if (edges_.element_list_sizes[i] > 2) {
      ++num_multi_updates;
    }
  }
  return num_multi_updates;
}

inline Int QuotientGraph::NumPivotCholeskyNonzeros() const QUOTIENT_NOEXCEPT {
  const Int pivot_size = std::abs(assembly_.signed_supernode_sizes[pivot_]);
  const Int structure_size =
      degree_lists_.degrees[pivot_] + assembly_.dense_supernode.size;
  const Int diag_block_nonzeros = (pivot_size * (pivot_size + 1)) / 2;
  const Int subdiagonal_nonzeros = structure_size * pivot_size;
  return diag_block_nonzeros + subdiagonal_nonzeros;
}

inline double QuotientGraph::NumPivotCholeskyFlops() const QUOTIENT_NOEXCEPT {
  const Int pivot_size = std::abs(assembly_.signed_supernode_sizes[pivot_]);
  const Int structure_size =
      degree_lists_.degrees[pivot_] + assembly_.dense_supernode.size;
  const double diag_block_flops = std::pow(1. * pivot_size, 3.) / 3.;
  const double schur_complement_flops =
      std::pow(1. * structure_size, 2.) * pivot_size;
  return diag_block_flops + schur_complement_flops;
}

inline void QuotientGraph::PackCountAndHashAdjacencies(
    Int i, Int num_elements, Int* degree, UInt* hash) QUOTIENT_NOEXCEPT {
  Int degree_new = *degree;
  UInt hash_new = *hash;

  Int num_packed = 0;
  const Int pack_index = edges_.element_list_offsets[i] + num_elements;
  const Int adjacency_list_beg =
      edges_.element_list_offsets[i] + edges_.element_list_sizes[i];
  QUOTIENT_ASSERT(pack_index <= adjacency_list_beg,
                  "Packing adjacencies after where they begin.");
  const Int adjacency_list_end =
      adjacency_list_beg + edges_.adjacency_list_sizes[i];
  for (Int k = adjacency_list_beg; k < adjacency_list_end; ++k) {
    const Int j = edges_.lists[k];
    // We must filter out: non-principal variables, eliminated supernodes,
    // the pivot supernode, and all members of the pivot structure.
    //
    // Due to our (temporary) negation of the supernode sizes of the pivot and
    // members of the pivot structure, all four cases can be handled by
    // demanding a positive (signed) supernode size.
    const Int supernode_size = assembly_.signed_supernode_sizes[j];
    if (supernode_size > 0) {
      degree_new += supernode_size;
      QUOTIENT_HASH_COMBINE(hash_new, j);
      edges_.lists[pack_index + num_packed++] = j;
    }
  }
  edges_.element_list_sizes[i] = num_elements;
  edges_.adjacency_list_sizes[i] = num_packed;
  *degree = degree_new;
  *hash = hash_new;
}

inline void QuotientGraph::InsertPivotElement(Int i) QUOTIENT_NOEXCEPT {
  Int& num_elements = edges_.element_list_sizes[i];
  const Int element_list_beg = edges_.element_list_offsets[i];
  const Int element_list_end = element_list_beg + num_elements;
  const Int adjacency_list_end =
      element_list_end + edges_.adjacency_list_sizes[i];
#ifdef QUOTIENT_DEBUG
  if (i == num_vertices_ - 1) {
    if (adjacency_list_end == Int(edges_.lists.Size())) {
      std::cerr << "Adjacency list ran off the end of the array." << std::endl;
    }
  } else {
    if (adjacency_list_end == edges_.element_list_offsets[i + 1]) {
      std::cerr << "Adjacency list overlapped with next element list."
                << std::endl;
    }
  }
#endif
  if (edges_.adjacency_list_sizes[i]) {
    // Move the first adjacency to the back.
    edges_.lists[adjacency_list_end] = edges_.lists[element_list_end];
  }
  edges_.lists[element_list_end] = pivot_;
  ++num_elements;
}

inline std::pair<Int, UInt> QuotientGraph::ExactEmptyDegreeAndHash(Int i)
    QUOTIENT_NOEXCEPT {
  const Int num_elements = edges_.element_list_sizes[i];
  Int degree = 0;
  UInt hash = 0;

  // We should only have one member of the element list, 'pivot'.
  QUOTIENT_ASSERT(
      num_elements == 0,
      "The element list was falsely assumed to have a single entry.");

  // Add |L_p \ supernode(i)|.
  const Int supernode_size = -assembly_.signed_supernode_sizes[i];
  QUOTIENT_ASSERT(supernode_size > 0,
                  "The flipped supernode size should have been positive.");
  degree += degree_lists_.degrees[pivot_] - supernode_size;
  QUOTIENT_HASH_COMBINE(hash, pivot_);

  // Add the cardinality of A_i \ supernode(i), where 'i' is the principal
  // variable.
  PackCountAndHashAdjacencies(i, num_elements, &degree, &hash);
  InsertPivotElement(i);

  return std::make_pair(degree, hash);
}

inline std::pair<Int, UInt> QuotientGraph::ExactSingleDegreeAndHash(Int i)
    QUOTIENT_NOEXCEPT {
  const Int offset = edges_.element_list_offsets[i];
  Int num_elements = edges_.element_list_sizes[i];
  Int degree = 0;
  UInt hash = 0;

  // There should be exactly two members in the element list, and the second
  // should be the pivot.
  QUOTIENT_ASSERT(num_elements == 1,
                  "The element list should have had exactly one members.");
  const Int& element = edges_.lists[offset];
  QUOTIENT_ASSERT(element != pivot_,
                  "The element list should have contained a non-pivot.");

  // Add |L_p \ supernode(i)|.
  const Int supernode_size = -assembly_.signed_supernode_sizes[i];
  QUOTIENT_ASSERT(supernode_size > 0,
                  "The flipped supernode size should have been positive.");
  degree += degree_lists_.degrees[pivot_] - supernode_size;
  QUOTIENT_HASH_COMBINE(hash, pivot_);

  Int& shifted_external_degree = node_flags_.flags[element];
  if (!shifted_external_degree) {
    --num_elements;
  } else {
    // Add |L_e \ L_p|.
    degree += shifted_external_degree - node_flags_.shift;
    QUOTIENT_HASH_COMBINE(hash, element);
  }

  // Add the cardinality of A_i \ supernode(i), where 'i' is the principal
  // variable.
  PackCountAndHashAdjacencies(i, num_elements, &degree, &hash);
  InsertPivotElement(i);

  return std::make_pair(degree, hash);
}

inline std::pair<Int, UInt> QuotientGraph::ExactGenericDegreeAndHash(Int i)
    QUOTIENT_NOEXCEPT {
  const Int offset = edges_.element_list_offsets[i];
  Int num_elements = edges_.element_list_sizes[i];
  Int degree = 0;
  UInt hash = 0;
  const Int shift = node_flags_.shift;

  // Add on the number of unique entries in the structures of the element lists
  // that are outside supernode(i).
  Int num_packed = 0;
  for (Int k = 0; k < num_elements; ++k) {
    const Int element = edges_.lists[offset + k];
    if (!node_flags_.flags[element]) {
      continue;
    }

    edges_.lists[offset + num_packed++] = element;
    QUOTIENT_HASH_COMBINE(hash, element);

    for (const Int& j : elements_[element]) {
      // Unabsorbed elements should not have any eliminated members of
      // in their element. Thus, we can take the absolute value of the signed
      // supernode size (since some members might be in the pivot structure).
      if (node_flags_.flags[j] == shift || i == j) {
        continue;
      }
      degree += std::abs(assembly_.signed_supernode_sizes[j]);
      node_flags_.flags[j] = shift;
    }
  }
  num_elements = num_packed;

  // Handle the pivot (but don't pack it yet).
  QUOTIENT_HASH_COMBINE(hash, pivot_);
  for (const Int& j : elements_[pivot_]) {
    // Unabsorbed elements should not have any eliminated members of
    // in their element. Thus, we can take the absolute value of the signed
    // supernode size (since some members might be in the pivot structure).
    if (node_flags_.flags[j] == shift || i == j) {
      continue;
    }
    degree += std::abs(assembly_.signed_supernode_sizes[j]);
    node_flags_.flags[j] = shift;
  }

  // Clear the mask.
  for (Int k = 0; k < num_elements; ++k) {
    const Int element = edges_.lists[offset + k];
    if (!node_flags_.flags[element]) {
      continue;
    }
    for (const Int& j : elements_[element]) {
      node_flags_.flags[j] = shift - 1;
    }
  }
  for (const Int& j : elements_[pivot_]) {
    node_flags_.flags[j] = shift - 1;
  }

  // Add the cardinality of A_i \ supernode(i), where 'i' is the principal
  // variable.
  PackCountAndHashAdjacencies(i, num_elements, &degree, &hash);
  InsertPivotElement(i);

  return std::make_pair(degree, hash);
}

inline void QuotientGraph::ExactDegreesAndHashes() QUOTIENT_NOEXCEPT {
  const Buffer<Int>& pivot_element = elements_[pivot_];
  const Int pivot_element_size = pivot_element.Size();

  for (Int index = 0; index < pivot_element_size; ++index) {
    const Int& i = pivot_element[index];

    std::pair<Int, UInt> degree_and_hash;
    const Int num_elements = edges_.element_list_sizes[i];
    if (num_elements == 0) {
      degree_and_hash = ExactEmptyDegreeAndHash(i);
    } else if (num_elements == 1) {
      degree_and_hash = ExactSingleDegreeAndHash(i);
    } else {
      degree_and_hash = ExactGenericDegreeAndHash(i);
    }
    degree_lists_.degrees[i] = degree_and_hash.first;
    hash_info_.lists.hashes[i] = degree_and_hash.second;
  }
}

inline void QuotientGraph::AmestoyDegreesAndHashes() QUOTIENT_NOEXCEPT {
  const Buffer<Int>& pivot_element = elements_[pivot_];
  const Int shift = node_flags_.shift;
  const Int pivot_degree = degree_lists_.degrees[pivot_];
  const Int pivot_element_size = pivot_element.Size();
  const Int num_vertices_left = num_vertices_ - num_eliminated_vertices_;

  for (Int index = 0; index < pivot_element_size; ++index) {
    const Int& i = pivot_element[index];
    Int degree = 0;
    UInt hash = 0;

    const Int supernode_size = -assembly_.signed_supernode_sizes[i];
    QUOTIENT_ASSERT(supernode_size > 0,
                    "The negated supernode size was expected to be positive.");
    const Int bound0 = num_vertices_left - supernode_size;

    // Note that this usage of 'external' refers to |L_p \ supernode(i)| and not
    // |L_e \ L_p|, as is the case for 'external_degrees'.
    const Int external_pivot_degree = pivot_degree - supernode_size;
    QUOTIENT_ASSERT(external_pivot_degree >= 0,
                    "Encountered a negative external pivot degree");

    const Int old_degree = degree_lists_.degrees[i];
    const Int bound1 = old_degree + external_pivot_degree;

    // bound_2 = |A_i \ supernode(i)| + |L_p \ supernode(i)| +
    //           \sum_{e in E_i \ {p}} |L_e \ L_p|.
    degree += external_pivot_degree;
    Int num_packed = 0;
    Int num_elements = edges_.element_list_sizes[i];
    const Int offset = edges_.element_list_offsets[i];
    for (Int k = offset; k < offset + num_elements; ++k) {
      const Int element = edges_.lists[k];
      QUOTIENT_ASSERT(element != pivot_, "Iterated over pivot element.");
      if (node_flags_.flags[element]) {
        QUOTIENT_ASSERT(
            node_flags_.flags[element],
            "Ran into an absorbed element in the external degree calculation.");
        edges_.lists[offset + num_packed++] = element;

        const Int external_degree = node_flags_.flags[element] - shift;
        QUOTIENT_ASSERT(external_degree >= 0, "Created a negative update.");
        degree += external_degree;
        QUOTIENT_HASH_COMBINE(hash, element);
      }
    }
    num_elements = num_packed;

    PackCountAndHashAdjacencies(i, num_elements, &degree, &hash);
    InsertPivotElement(i);
    QUOTIENT_HASH_COMBINE(hash, pivot_);

    degree = std::min(degree, bound0);
    degree = std::min(degree, bound1);

    degree_lists_.degrees[i] = degree;
    hash_info_.lists.hashes[i] = hash;
  }
}

inline void QuotientGraph::AshcraftDegreesAndHashes() QUOTIENT_NOEXCEPT {
  const Buffer<Int>& pivot_element = elements_[pivot_];
  const Int pivot_element_size = pivot_element.Size();

  for (Int index = 0; index < pivot_element_size; ++index) {
    const Int& i = pivot_element[index];

    std::pair<Int, UInt> degree_and_hash;
    // Note that the list size is one *before* adding the pivot.
    if (edges_.element_list_sizes[i] == 1) {
      degree_and_hash = ExactSingleDegreeAndHash(i);
    } else {
      degree_and_hash = GilbertDegreeAndHash(i);
    }
    degree_lists_.degrees[i] = degree_and_hash.first;
    hash_info_.lists.hashes[i] = degree_and_hash.second;
  }
}

inline std::pair<Int, UInt> QuotientGraph::GilbertDegreeAndHash(Int i)
    QUOTIENT_NOEXCEPT {
  Int degree = 0;
  UInt hash = 0;

  const Int supernode_size = -assembly_.signed_supernode_sizes[i];
  QUOTIENT_ASSERT(supernode_size > 0,
                  "The negated supernode size was expected to be positive.");

  Int num_packed = 0;
  const Int offset = edges_.element_list_offsets[i];
  Int num_elements = edges_.element_list_sizes[i];
  for (Int k = 0; k < num_elements; ++k) {
    const Int element = edges_.lists[offset + k];
    if (!node_flags_.flags[element]) {
      continue;
    }
    edges_.lists[offset + num_packed++] = element;

    QUOTIENT_ASSERT(assembly_.parent_or_tail[element] < 0,
                    "Used absorbed element in Gilbert degree update.");
    QUOTIENT_ASSERT(degree_lists_.degrees[element] - supernode_size >= 0,
                    "Negative Gilbert degree update.");
    degree += degree_lists_.degrees[element] - supernode_size;
    QUOTIENT_HASH_COMBINE(hash, element);
  }
  num_elements = num_packed;
  degree += degree_lists_.degrees[pivot_] - supernode_size;
  QUOTIENT_HASH_COMBINE(hash, pivot_);

  PackCountAndHashAdjacencies(i, num_elements, &degree, &hash);
  InsertPivotElement(i);

  const Int num_vertices_left = num_vertices_ - num_eliminated_vertices_;
  const Int bound = num_vertices_left - supernode_size;
  degree = std::min(degree, bound);

  return std::make_pair(degree, hash);
}

inline void QuotientGraph::GilbertDegreesAndHashes() QUOTIENT_NOEXCEPT {
  const Buffer<Int>& pivot_element = elements_[pivot_];
  const Int num_vertices_left = num_vertices_ - num_eliminated_vertices_;
  const Int pivot_degree = degree_lists_.degrees[pivot_];
  const Int pivot_element_size = pivot_element.Size();

  for (Int index = 0; index < pivot_element_size; ++index) {
    const Int& i = pivot_element[index];
    Int degree = 0;
    UInt hash = 0;

    const Int supernode_size = -assembly_.signed_supernode_sizes[i];
    QUOTIENT_ASSERT(supernode_size > 0,
                    "The negated supernode size was expected to be positive.");

    Int num_packed = 0;
    Int num_elements = edges_.element_list_sizes[i];
    const Int offset = edges_.element_list_offsets[i];
    for (Int k = 0; k < num_elements; ++k) {
      const Int element = edges_.lists[offset + k];
      if (!node_flags_.flags[element]) {
        continue;
      }
      edges_.lists[offset + num_packed++] = element;

      QUOTIENT_ASSERT(assembly_.parent_or_tail[element] < 0,
                      "Used absorbed element in Gilbert degree update.");
      QUOTIENT_ASSERT(degree_lists_.degrees[element] - supernode_size >= 0,
                      "Negative Gilbert degree update.");
      degree += degree_lists_.degrees[element] - supernode_size;
      QUOTIENT_HASH_COMBINE(hash, element);
    }
    num_elements = num_packed;
    degree += pivot_degree - supernode_size;
    QUOTIENT_HASH_COMBINE(hash, pivot_);

    PackCountAndHashAdjacencies(i, num_elements, &degree, &hash);
    InsertPivotElement(i);

    const Int bound = num_vertices_left - supernode_size;
    degree = std::min(degree, bound);

    degree_lists_.degrees[i] = degree;
    hash_info_.lists.hashes[i] = hash;
  }
}

inline void QuotientGraph::ComputeDegreesAndHashes() QUOTIENT_NOEXCEPT {
  QUOTIENT_START_TIMER(timers_, kComputeDegrees);

  // Remove the old values from the linked lists.
  for (const Int& i : elements_[pivot_]) {
    degree_lists_.RemoveDegree(i);
  }

  switch (control_.degree_type) {
    case kExactDegree: {
      ExactDegreesAndHashes();
      break;
    }
    case kAmestoyDegree: {
      AmestoyDegreesAndHashes();
      break;
    }
    case kAshcraftDegree: {
      AshcraftDegreesAndHashes();
      break;
    }
    case kGilbertDegree: {
      GilbertDegreesAndHashes();
      break;
    }
  }

  // Modify the linked lists to contain the new values.
  for (const Int& i : elements_[pivot_]) {
    const Int degree = degree_lists_.degrees[i];
    degree_lists_.AddDegree(i, degree);
  }

  QUOTIENT_STOP_TIMER(timers_, kComputeDegrees);
}

inline bool QuotientGraph::StructuralVariablesAreQuotientIndistinguishable(
    Int i, Int j) const QUOTIENT_NOEXCEPT {
  QUOTIENT_ASSERT(
      i != j,
      "Explicitly tested for equivalence of a supernode with itself...");
  const Int offset_i = edges_.element_list_offsets[i];
  const Int offset_j = edges_.element_list_offsets[j];
  const Int num_elements_i = edges_.element_list_sizes[i];
  const Int num_elements_j = edges_.element_list_sizes[j];
  const Int num_adjacencies_i = edges_.adjacency_list_sizes[i];
  const Int num_adjacencies_j = edges_.adjacency_list_sizes[j];

#ifdef QUOTIENT_DEBUG
  for (Int k = 0; k < num_elements_i; ++k) {
    const Int element = edges_.lists[offset_i + k];
    QUOTIENT_ASSERT(
        assembly_.parent_or_tail[element] < 0,
        "Absorbed element was in element list during absorption test.");
  }
  for (Int k = 0; k < num_elements_j; ++k) {
    const Int element = edges_.lists[offset_j + k];
    QUOTIENT_ASSERT(
        assembly_.parent_or_tail[element] < 0,
        "Absorbed element was in element list during absorption test.");
  }
  for (Int k = 0; k < num_adjacencies_i; ++k) {
    const Int index = edges_.lists[offset_i + num_elements_i + k];
    QUOTIENT_ASSERT(
        assembly_.signed_supernode_sizes[index] > 0,
        "Non-positive supernode size in adj list during absorption test.");
  }
  for (Int k = 0; k < num_adjacencies_j; ++k) {
    const Int index = edges_.lists[offset_j + num_elements_j + k];
    QUOTIENT_ASSERT(
        assembly_.signed_supernode_sizes[index] > 0,
        "Non-positive supernode size in adj list during absorption test.");
  }
#endif

  // We must have the same number of elements and adjacencies to have a match.
  if (num_elements_i != num_elements_j ||
      num_adjacencies_i != num_adjacencies_j) {
    return false;
  }

  // Check if E_i = E_j.
  for (Int k = 0; k < num_elements_i; ++k) {
    if (edges_.lists[offset_i + k] != edges_.lists[offset_j + k]) {
      return false;
    }
  }

  // Check if A_i = A_j by ensuring that all members of A_j are in A_i (by
  // checking that the indices are currently flagged).
  const Int shift = node_flags_.shift;
  const Int adjacency_list_j_beg = offset_j + num_elements_j;
  const Int adjacency_list_j_end = adjacency_list_j_beg + num_adjacencies_j;
  for (Int k = adjacency_list_j_beg; k < adjacency_list_j_end; ++k) {
    const Int index = edges_.lists[k];
    if (node_flags_.flags[index] != shift) {
      return false;
    }
  }

  return true;
}

inline void QuotientGraph::MergeVariables() QUOTIENT_NOEXCEPT {
  QUOTIENT_START_TIMER(timers_, kMergeVariables);
  const Buffer<Int>& pivot_element = elements_[pivot_];
  const Int supernodal_struct_size = pivot_element.Size();
  const Int shift = node_flags_.shift;

  // Add the hashes into the hash lists.
  for (Int i_index = 0; i_index < supernodal_struct_size; ++i_index) {
    const Int i = pivot_element[i_index];
    const UInt hash = hash_info_.lists.hashes[i];
    const Int bucket = hash % num_vertices_;
    hash_info_.lists.AddHash(i, hash, bucket);
  }

  Int num_merges = 0;
  const Buffer<Int>& next_member = hash_info_.lists.next_member;
  for (Int i_index = 0; i_index < supernodal_struct_size; ++i_index) {
    Int i = pivot_element[i_index];
    const Int bucket = hash_info_.lists.buckets[i];
    const Int head = hash_info_.lists.heads[bucket];
    // We are not the head of the bucket.
    if (i != head) {
      continue;
    }

    // Test the unique pairs in the bucket.
    for (; next_member[i] != -1; i = next_member[i]) {
      if (!assembly_.signed_supernode_sizes[i]) {
        continue;
      }
      QUOTIENT_ASSERT(assembly_.signed_supernode_sizes[i] < 0,
                      "Supernode size should have been temporarily negative.");
      const UInt i_hash = hash_info_.lists.hashes[i];
      bool scattered_adjacencies = false;
      const Int adjacency_list_i_beg =
          edges_.element_list_offsets[i] + edges_.element_list_sizes[i];
      const Int adjacency_list_i_end =
          adjacency_list_i_beg + edges_.adjacency_list_sizes[i];
      for (Int j = next_member[i]; j != -1; j = next_member[j]) {
        if (!assembly_.signed_supernode_sizes[j]) {
          continue;
        }
        QUOTIENT_ASSERT(
            assembly_.signed_supernode_sizes[j] < 0,
            "Supernode size should have been temporarily negative.");
        const UInt j_hash = hash_info_.lists.hashes[j];
        if (i_hash != j_hash) {
          ++hash_info_.num_bucket_collisions;
          continue;
        }
        if (!scattered_adjacencies) {
          // Temporarily flag the adjacencies of supervariable i by setting
          // 'node_flags_.flags' to the current shift value in said indices.
          // This is representative because only elements should hold this
          // value.
          //
          // It's worth noting that the adjacency lists are not guaranteed to
          // be in the same order (the reason for the scattering) because of the
          // mechanism used to insert the pivot element into the element list:
          // it shifts the first adjacency to the back of the list. Different
          // variables could have done so in a different manner and (through
          // adjacency removals) ended up with the same adjacency set in
          // different orders.
          //
          // The element lists are unsorted but should always appear in the
          // same order.
          //
          // TODO(Jack Poulson): Consider whether it would make more sense to
          // sort the adjacency lists and then do the trivial comparison. This
          // would have a significant advantage if multithreaded parallelism
          // were to be introduced into MergeVariables, as it would eliminate
          // the need for each thread to have a separate 'flags' buffer.
          for (Int k = adjacency_list_i_beg; k < adjacency_list_i_end; ++k) {
            const Int index = edges_.lists[k];
            node_flags_.flags[index] = shift;
          }
          scattered_adjacencies = true;
        }

        if (StructuralVariablesAreQuotientIndistinguishable(i, j)) {
          ++num_merges;
          const Int absorbed_size = -assembly_.signed_supernode_sizes[j];
          QUOTIENT_ASSERT(absorbed_size > 0,
                          "Absorbed size should have been positive.");

          const Int old_degree = degree_lists_.degrees[i];
          const Int degree = old_degree - absorbed_size;
          degree_lists_.UpdateDegree(i, degree);
          degree_lists_.RemoveDegree(j);

          // Merge [i] -> [j] (and i becomes the principal member).
          //
          // Recall that the supernode sizes of principal variables in the
          // pivot structure have been negated.
          const Int tail_index_j = SYMMETRIC_INDEX(assembly_.parent_or_tail[j]);
          assembly_.parent_or_tail[j] = i;
          assembly_.parent_or_tail[i] = SYMMETRIC_INDEX(tail_index_j);
          assembly_.signed_supernode_sizes[i] -= absorbed_size;
          assembly_.signed_supernode_sizes[j] = 0;

          edges_.element_list_sizes[j] = 0;
          edges_.adjacency_list_sizes[j] = 0;
        } else {
          ++hash_info_.num_collisions;
          ++hash_info_.num_bucket_collisions;
        }
      }
      if (scattered_adjacencies) {
        // Reset the flagged adjacencies to the shift minus one.
        for (Int k = adjacency_list_i_beg; k < adjacency_list_i_end; ++k) {
          const Int index = edges_.lists[k];
          node_flags_.flags[index] = shift - 1;
        }
      }
    }
    hash_info_.lists.ClearBucket(bucket);
  }

#ifdef QUOTIENT_DEBUG
  for (Int index = 0; index < num_vertices_; ++index) {
    if (hash_info_.lists.heads[index] != -1) {
      std::cerr << "Did not clear head for bucket " << index << std::endl;
    }
  }
#endif

  QUOTIENT_STOP_TIMER(timers_, kMergeVariables);
}

inline void QuotientGraph::FinalizePivot() QUOTIENT_NOEXCEPT {
  QUOTIENT_START_TIMER(timers_, kFinalizePivot);
  ResetExternalDegrees();

  const Int supernode_size = -assembly_.signed_supernode_sizes[pivot_];
  QUOTIENT_ASSERT(supernode_size > 0,
                  "The supernode size was assumed positive.");

  // Since this supervariable is being eliminated, it needs to be implicitly
  // removed from all elements containing it. This is accomplished through
  // decrementing their element sizes and marking the supernode as eliminated
  // by making its sign negative.
  const Int offset = edges_.element_list_offsets[pivot_];
  const Int num_elements = edges_.element_list_sizes[pivot_];
  for (Int k = offset; k < offset + num_elements; ++k) {
    const Int element = edges_.lists[k];
    degree_lists_.degrees[element] -= supernode_size;
  }

  // Unflip the signs of supernode sizes for the members of the pivot structure
  // and simultaneously remove the non-principal members.
  Int num_packed = 0;
  for (const Int& i : elements_[pivot_]) {
    QUOTIENT_ASSERT(
        assembly_.signed_supernode_sizes[i] <= 0,
        "A member of pivot element had a positive signed supernode size.");
    if (assembly_.signed_supernode_sizes[i] < 0) {
      assembly_.signed_supernode_sizes[i] =
          -assembly_.signed_supernode_sizes[i];
      elements_[pivot_][num_packed++] = i;
    }
  }
  elements_[pivot_].Resize(num_packed);

  elimination_order_.push_back(pivot_);
  num_eliminated_vertices_ += supernode_size;

  QUOTIENT_STOP_TIMER(timers_, kFinalizePivot);
}

inline void QuotientGraph::ComputePostorder(Buffer<Int>* postorder) const
    QUOTIENT_NOEXCEPT {
  auto supernode_principal = [&](Int index) {
    while (!assembly_.signed_supernode_sizes[index]) {
      QUOTIENT_ASSERT(
          assembly_.parent_or_tail[index] >= 0,
          "Negative member of assembly_.parent_or_tail while computing "
          "supernode principal.");
      index = assembly_.parent_or_tail[index];
    }
    return index;
  };

  // Fill the supernode non-principal member lists.
  Buffer<Int> nonprincipal_offsets(num_vertices_ + 1, 0);
  Int num_nonprincipal_members = 0;
  for (Int i = 0; i < num_vertices_; ++i) {
    nonprincipal_offsets[i] = num_nonprincipal_members;
    const Int supernode_size = -assembly_.signed_supernode_sizes[i];
    QUOTIENT_ASSERT(supernode_size >= 0, "Supernode size was negative.");
    if (supernode_size) {
      num_nonprincipal_members += supernode_size - 1;
    }
  }
  nonprincipal_offsets[num_vertices_] = num_nonprincipal_members;
  auto offsets_copy = nonprincipal_offsets;
  Buffer<Int> nonprincipal_members(num_nonprincipal_members);
  for (Int i = 0; i < num_vertices_; ++i) {
    const Int supernode_size = -assembly_.signed_supernode_sizes[i];
    if (!supernode_size) {
      const Int principal = supernode_principal(i);
      nonprincipal_members[offsets_copy[principal]++] = i;
    }
  }

  // Reconstruct the child links from the parent links in a contiguous array
  // (similar to the CSR format) by first counting the number of children of
  // each node.
  Buffer<Int> children;
  Buffer<Int> child_offsets;
  ChildrenFromParentSubsequence(assembly_.parent_or_tail, elimination_order_,
                                &children, &child_offsets);

  // Scan for the roots and launch a pe-order traversal on each of them.
  // We march through elimination_order in reverse order so that, after a
  // subsequent reversal, the lowest degree nodes are in the upper-left.
  postorder->Resize(num_vertices_);
  Int* iter = postorder->begin();
  for (const Int& i : elimination_order_) {
    if (assembly_.parent_or_tail[i] >= 0) {
      // This element was absorbed into another element.
      continue;
    }
    iter = PreorderTree(i, nonprincipal_members, nonprincipal_offsets, children,
                        child_offsets, iter);
  }

  // Reverse the preordering (to form a postordering) in-place.
  std::reverse(postorder->begin(), iter);

  if (assembly_.dense_supernode.size) {
    // Pack the dense columns into the back.
    for (Int i = 0; i < num_vertices_; ++i) {
      if (!assembly_.signed_supernode_sizes[i] &&
          assembly_.parent_or_tail[i] < 0) {
        *iter = i;
        ++iter;
      }
    }
  }

  QUOTIENT_ASSERT(iter == postorder->end(),
                  "Postorder had incorrect final offset.");
}

inline Int* QuotientGraph::PreorderTree(Int root,
                                        const Buffer<Int>& nonprincipal_members,
                                        const Buffer<Int>& nonprincipal_offsets,
                                        const Buffer<Int>& children,
                                        const Buffer<Int>& child_offsets,
                                        Int* iter) const QUOTIENT_NOEXCEPT {
  std::vector<Int> stack;
  stack.reserve(num_vertices_);

  stack.push_back(root);

  while (!stack.empty()) {
    // Pop a principal variable from the stack.
    const Int element = stack.back();
    stack.pop_back();

    // Push the supernode into the preorder in reverse order so that, when we
    // later call std::reverse to generate a postorder, the principal member of
    // the supernode comes first.
    const Int nonprincipal_beg = nonprincipal_offsets[element];
    const Int nonprincipal_end = nonprincipal_offsets[element + 1];
    for (Int j = nonprincipal_end - 1; j >= nonprincipal_beg; --j) {
      (*iter++) = nonprincipal_members[j];
    }
    *(iter++) = element;

    // Push the children onto the stack.
    const Int child_beg = child_offsets[element];
    const Int child_end = child_offsets[element + 1];
    for (Int index = child_beg; index < child_end; ++index) {
      stack.push_back(children[index]);
    }
  }

  return iter;
}

inline void QuotientGraph::PermutedSupernodeSizes(
    const Buffer<Int>& inverse_permutation,
    Buffer<Int>* permuted_supernode_sizes) const QUOTIENT_NOEXCEPT {
  const Int num_supernodes = elimination_order_.size();
  permuted_supernode_sizes->Clear();
  permuted_supernode_sizes->Resize(num_supernodes);

  Int i_perm = 0;
  for (Int index = 0; index < num_supernodes; ++index) {
    const Int i = inverse_permutation[i_perm];
    QUOTIENT_ASSERT(assembly_.signed_supernode_sizes[i] < 0,
                    "Supernode size was negative.");
    const Int supernode_size = -assembly_.signed_supernode_sizes[i];
    (*permuted_supernode_sizes)[index] = supernode_size;

    i_perm += supernode_size;
  }
}

inline void QuotientGraph::PermutedMemberToSupernode(
    const Buffer<Int>& inverse_permutation,
    Buffer<Int>* permuted_member_to_supernode) const QUOTIENT_NOEXCEPT {
  const Int num_indices = inverse_permutation.Size();
  permuted_member_to_supernode->Clear();
  permuted_member_to_supernode->Resize(num_indices);

  Int permuted_supernode = -1;
  for (Int i_perm = 0; i_perm < num_indices; ++i_perm) {
    const Int i = inverse_permutation[i_perm];
    if (assembly_.signed_supernode_sizes[i]) {
      // 'i' is a principal variable of an absorbed element.
      QUOTIENT_ASSERT(assembly_.signed_supernode_sizes[i] < 0,
                      "Supernode size was negative.");
      (*permuted_member_to_supernode)[i_perm] = ++permuted_supernode;
    } else {
      QUOTIENT_ASSERT(permuted_supernode >= 0, "Invalid permuted ordering.");
      (*permuted_member_to_supernode)[i_perm] = permuted_supernode;
    }
  }
}

inline void QuotientGraph::PermutedAssemblyParents(
    const Buffer<Int>& permutation,
    const Buffer<Int>& permuted_member_to_supernode,
    Buffer<Int>* permuted_assembly_parents) const QUOTIENT_NOEXCEPT {
  const Int num_supernodes = elimination_order_.size();
  permuted_assembly_parents->Clear();
  permuted_assembly_parents->Resize(num_supernodes);
  for (Int index = 0; index < num_supernodes; ++index) {
    const Int original_principal = elimination_order_[index];
    const Int original_parent = assembly_.parent_or_tail[original_principal];

    const Int permuted_principal = permutation[original_principal];
    const Int permuted_parent =
        original_parent >= 0 ? permutation[original_parent] : -1;

    const Int permuted_principal_supernode =
        permuted_member_to_supernode[permuted_principal];
    const Int permuted_parent_supernode =
        permuted_parent >= 0 ? permuted_member_to_supernode[permuted_parent]
                             : -1;

    (*permuted_assembly_parents)[permuted_principal_supernode] =
        permuted_parent_supernode;
  }
}

inline void QuotientGraph::ExternalDegrees() QUOTIENT_NOEXCEPT {
  QUOTIENT_START_TIMER(timers_, kExternalDegrees);
  const Int shift = node_flags_.shift;
  const bool aggressive_absorption = control_.aggressive_absorption;

  for (const Int& i : elements_[pivot_]) {
    const Int supernode_size = -assembly_.signed_supernode_sizes[i];
    QUOTIENT_ASSERT(supernode_size > 0,
                    "supernode " + std::to_string(i) +
                        " had non-positive signed size "
                        "when computing external element sizes");
    const Int shift_minus_supernode_size = shift - supernode_size;

    const Int offset = edges_.element_list_offsets[i];
    const Int num_elements = edges_.element_list_sizes[i];
    for (Int k = offset; k < offset + num_elements; ++k) {
      const Int element = edges_.lists[k];
      Int& shifted_external_degree = node_flags_.flags[element];

      if (shifted_external_degree >= shift) {
        QUOTIENT_ASSERT(assembly_.parent_or_tail[element] < 0,
                        "Tried to subtract from an absorbed element.");
        QUOTIENT_ASSERT(
            shifted_external_degree != shift,
            "Shifted external size was equal to shift before subtracting");
        shifted_external_degree -= supernode_size;
      } else if (shifted_external_degree) {
        QUOTIENT_ASSERT(assembly_.parent_or_tail[element] < 0,
                        "Tried to subtract from an absorbed element.");
        shifted_external_degree =
            degree_lists_.degrees[element] + shift_minus_supernode_size;
      }
      QUOTIENT_ASSERT(
          !shifted_external_degree || shifted_external_degree >= shift,
          "Computed negative external element degree.");
      if (aggressive_absorption && shifted_external_degree == shift) {
        ++num_aggressive_absorptions_;
        shifted_external_degree = 0;
        assembly_.parent_or_tail[element] = pivot_;
        elements_[element].Clear();
        edges_.element_list_sizes[element] = 0;
      }
    }
  }

  QUOTIENT_STOP_TIMER(timers_, kExternalDegrees);
}

inline void QuotientGraph::ResetExternalDegrees() QUOTIENT_NOEXCEPT {
  if (node_flags_.shift + node_flags_.max_degree < node_flags_.shift_cap) {
    node_flags_.shift += node_flags_.max_degree + 1;
    return;
  }

// Reset all non-absorbed shifted degrees to 1. Afterwards, the absorbed
// elements will be marked as 0, the unabsorbed as 1, and the shift will be
// 2.
#ifdef QUOTIENT_DEBUG
  std::cerr << "Resetting external degrees." << std::endl;
#endif
  node_flags_.shift = 2;
  for (Int i = 0; i < node_flags_.flags.Size(); ++i) {
    if (node_flags_.flags[i]) {
      node_flags_.flags[i] = 1;
    }
  }
}

inline void QuotientGraph::CombineDenseNodes() QUOTIENT_NOEXCEPT {
  if (!assembly_.dense_supernode.size) {
    return;
  }

  // Absorb all of the dense nodes into the principal member.
  for (Int i = 0; i < num_vertices_; ++i) {
    if (!assembly_.signed_supernode_sizes[i] &&
        assembly_.parent_or_tail[i] < 0) {
      if (i == assembly_.dense_supernode.principal_member) {
        assembly_.signed_supernode_sizes[i] = -assembly_.dense_supernode.size;
      } else {
        assembly_.parent_or_tail[i] =
            assembly_.dense_supernode.principal_member;
        assembly_.signed_supernode_sizes[i] = 0;
      }
    }
  }

  // Point the non-dense elements currently marked as roots to the dense
  // supernode.
  for (Int i = 0; i < num_vertices_; ++i) {
    if (assembly_.signed_supernode_sizes[i] &&
        i != assembly_.dense_supernode.principal_member &&
        assembly_.parent_or_tail[i] < 0) {
      assembly_.parent_or_tail[i] = assembly_.dense_supernode.principal_member;
    }
  }
}

inline Buffer<std::pair<std::string, double>> QuotientGraph::ComponentSeconds()
    const QUOTIENT_NOEXCEPT {
  Buffer<std::pair<std::string, double>> times;
#ifdef QUOTIENT_ENABLE_TIMERS
  for (const std::pair<std::string, Timer>& pairing : timers_) {
    times.emplace_back(pairing.first, pairing.second.TotalSeconds());
  }
#endif
  return times;
}

}  // namespace quotient

#endif  // ifndef QUOTIENT_QUOTIENT_GRAPH_IMPL_H_
