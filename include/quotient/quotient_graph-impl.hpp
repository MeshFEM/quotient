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

static constexpr char kSetup[] = "Setup";
static constexpr char kComputePivotStructure[] = "ComputePivotStructure";
static constexpr char kExternalDegrees[] = "ExternalDegrees";
static constexpr char kComputeDegreesAndHashes[] = "ComputeDegreesAndHashes";
static constexpr char kMergeVariables[] = "MergeVariables";
static constexpr char kFinalizePivot[] = "FinalizePivot";
static constexpr char kComputePostorder[] = "ComputePostorder";
#else
#define QUOTIENT_START_TIMER(timer, name)
#define QUOTIENT_STOP_TIMER(timer, name)
#endif  // ifdef QUOTIENT_ENABLE_TIMERS

inline bool QuotientGraph::QuotientGraphData::ActiveSupernode(Int i) const
    QUOTIENT_NOEXCEPT {
  return element_offsets[i] >= -1;
}

inline Int QuotientGraph::QuotientGraphData::Parent(Int i) const
    QUOTIENT_NOEXCEPT {
  QUOTIENT_ASSERT(element_offsets[i] < 0,
                  "Cannot retrieve parent of active object.");
  return SYMMETRIC_INDEX(element_offsets[i]);
}

inline void QuotientGraph::QuotientGraphData::SetParent(Int i, Int parent)
    QUOTIENT_NOEXCEPT {
  element_offsets[i] = SYMMETRIC_INDEX(parent);
}

inline const Int* QuotientGraph::QuotientGraphData::ElementList(Int i) const
    QUOTIENT_NOEXCEPT {
  return &lists[element_offsets[i]];
}

inline Int* QuotientGraph::QuotientGraphData::ElementList(Int i)
    QUOTIENT_NOEXCEPT {
  return &lists[element_offsets[i]];
}

inline const Int* QuotientGraph::QuotientGraphData::ElementData(Int i) const
    QUOTIENT_NOEXCEPT {
  return &lists[element_offsets[i]];
}

inline Int* QuotientGraph::QuotientGraphData::ElementData(Int i)
    QUOTIENT_NOEXCEPT {
  return &lists[element_offsets[i]];
}

inline const Int* QuotientGraph::QuotientGraphData::AdjacencyList(Int i) const
    QUOTIENT_NOEXCEPT {
  return &lists[element_offsets[i] + element_sizes[i]];
}

inline Int* QuotientGraph::QuotientGraphData::AdjacencyList(Int i)
    QUOTIENT_NOEXCEPT {
  return &lists[element_offsets[i] + element_sizes[i]];
}

inline const Int& QuotientGraph::QuotientGraphData::ElementListSize(
    Int element) const QUOTIENT_NOEXCEPT {
  return element_sizes[element];
}

inline Int& QuotientGraph::QuotientGraphData::ElementListSize(Int element)
    QUOTIENT_NOEXCEPT {
  return element_sizes[element];
}

inline const Int& QuotientGraph::QuotientGraphData::ElementSize(
    Int element) const QUOTIENT_NOEXCEPT {
  return element_sizes[element];
}

inline Int& QuotientGraph::QuotientGraphData::ElementSize(Int element)
    QUOTIENT_NOEXCEPT {
  return element_sizes[element];
}

inline void QuotientGraph::QuotientGraphData::Pack() QUOTIENT_NOEXCEPT {
  const Int offset_save = offset;

  // Overwrite the offsets of all of the active variables and elements
  // with the first entry of the object and replace the entry with the
  // symmetric version of the object index (as a flag).
  const Int num_vertices = element_sizes.Size();
  for (Int i = 0; i < num_vertices; ++i) {
    if (ActiveSupernode(i)) {
      QUOTIENT_ASSERT(element_sizes[i] + adjacency_list_sizes[i] != 0,
                      "Had a zero-length active variable.");
      Int* element_list = &lists[element_offsets[i]];
      element_offsets[i] = element_list[0];
      element_list[0] = SYMMETRIC_INDEX(i);
    }
  }

  // Pack the adjacencies and elements.
  Int pack_offset = 0;
  Int read_offset = 0;
  while (read_offset < offset_save) {
    const Int entry = lists[read_offset++];
    if (entry >= 0) {
      continue;
    }
    const Int i = SYMMETRIC_INDEX(entry);
    const Int element_size = element_sizes[i];
    const Int adjacency_size = adjacency_list_sizes[i];
    QUOTIENT_ASSERT(element_size || adjacency_size, "Packing empty element");

    // The current index is the beginning of an object.
    const Int displaced_entry = element_offsets[i];
    element_offsets[i] = pack_offset;
    lists[pack_offset++] = displaced_entry;

    // Pack the rest of the object.
    const Int length = element_size + adjacency_size;
    for (Int k = 1; k < length; ++k) {
      lists[pack_offset++] = lists[read_offset++];
    }
  }

  offset = pack_offset;
#ifdef QUOTIENT_DEBUG
  std::cout << "Began with " << offset_save << ", ended with " << pack_offset
            << std::endl;
#endif
}

inline QuotientGraph::QuotientGraph(const CoordinateGraph& graph,
                                    const MinimumDegreeControl& control)
    : QuotientGraph(graph.NumSources(), graph.Edges(), control) {}

inline Int QuotientGraph::ConvertEdgeCountsIntoOffsets() QUOTIENT_NOEXCEPT {
  // Every row with at least this many non-diagonal nonzeros will be treated
  // as dense and moved to the back of the postordering.
  const Int dense_threshold =
      std::max(1.f * control_.min_dense_threshold,
               control_.dense_sqrt_multiple * std::sqrt(1.f * num_vertices_));

  // Convert the counts into an offset scan.
  Int num_edges = 0;
  for (Int source = 0; source < num_vertices_; ++source) {
    graph_data_.element_offsets[source] = num_edges;
    if (graph_data_.adjacency_list_sizes[source] >= dense_threshold) {
      if (!graph_data_.dense_supernode.size) {
        graph_data_.dense_supernode.principal_member = source;
        ++num_eliminated_supernodes_;
      }
      ++graph_data_.dense_supernode.size;
      ++num_eliminated_vertices_;
      // We will denote this source as dense by setting its supernode size to
      // 0 and keeping its parent as -1. A merged variable also has a supernode
      // size of 0, but its parent is a valid vertex index.
      graph_data_.signed_supernode_sizes[source] = 0;
      graph_data_.adjacency_list_sizes[source] = 0;
    }
    num_edges += graph_data_.adjacency_list_sizes[source];
  }
#ifdef QUOTIENT_DEBUG
  if (graph_data_.dense_supernode.size) {
    std::cout << "Eliminated " << graph_data_.dense_supernode.size
              << " dense rows." << std::endl;
  }
#endif

  return num_edges;
}

inline void QuotientGraph::InitializeDegreeAndHashLists() QUOTIENT_NOEXCEPT {
  degrees_and_hashes_.lists.degrees.Resize(num_vertices_, 0);
  degrees_and_hashes_.lists.heads.Resize(num_vertices_, -1);
  degrees_and_hashes_.lists.next_member.Resize(num_vertices_, -1);
  degrees_and_hashes_.lists.last_member.Resize(num_vertices_, -1);
  for (Int source = 0; source < num_vertices_; ++source) {
    if (!graph_data_.signed_supernode_sizes[source]) {
      // Skip the dense row.
      continue;
    }
    const Int degree = graph_data_.adjacency_list_sizes[source];
    degrees_and_hashes_.lists.AddDegree(source, degree);
  }
  degrees_and_hashes_.num_collisions = 0;
  degrees_and_hashes_.num_bucket_collisions = 0;
}

inline void QuotientGraph::InitializeNodeFlags() QUOTIENT_NOEXCEPT {
  // The absorbed elements will be maintained as 0, the unabsorbed will be
  // initialized as 1, and the shift will be initialized as 2.
  node_flags_.shift = 2;
  node_flags_.flags.Resize(num_vertices_, 1);
  node_flags_.max_degree = 0;
  node_flags_.shift_cap = std::numeric_limits<Int>::max() - num_vertices_;
}

inline QuotientGraph::QuotientGraph(Int num_vertices,
                                    const Buffer<GraphEdge>& edges,
                                    const MinimumDegreeControl& control)
    : control_(control),
      num_vertices_(num_vertices),
      num_eliminated_vertices_(0),
      num_eliminated_supernodes_(0),
      num_aggressive_absorptions_(0) {
  QUOTIENT_START_TIMER(timers_, kSetup);

  graph_data_.signed_supernode_sizes.Resize(num_vertices_, 1);
  graph_data_.dense_supernode.size = 0;
  graph_data_.dense_supernode.principal_member = -1;

  // Initialize the counts for the adjacency lists.
  graph_data_.element_offsets.Resize(num_vertices_);
  graph_data_.element_sizes.Resize(num_vertices_, 0);
  graph_data_.adjacency_list_sizes.Resize(num_vertices_, 0);
  for (const GraphEdge& edge : edges) {
    if (edge.first != edge.second) {
      ++graph_data_.adjacency_list_sizes[edge.first];
    }
  }

  Int num_edges = ConvertEdgeCountsIntoOffsets();

  // Pack the edges. We leave extra room for packing element structures.
  // TODO(Jack Poulson): Make this coefficient configurable.
  const float kNumEdgesRatio = 0.2f;
  const Int extra_element_space = kNumEdgesRatio * num_edges + num_vertices_;
  graph_data_.lists.Resize(num_edges + extra_element_space);
  num_edges = 0;
  for (const GraphEdge& edge : edges) {
    const Int source = edge.first;
    if (graph_data_.signed_supernode_sizes[source] && edge.second != source) {
      graph_data_.lists[num_edges++] = edge.second;
    }
  }
  graph_data_.offset = num_edges;

  InitializeDegreeAndHashLists();
  InitializeNodeFlags();

  QUOTIENT_STOP_TIMER(timers_, kSetup);
}

template <typename Field>
inline QuotientGraph::QuotientGraph(Int num_vertices,
                                    const Buffer<MatrixEntry<Field>>& entries,
                                    const MinimumDegreeControl& control)
    : control_(control),
      num_vertices_(num_vertices),
      num_eliminated_vertices_(0),
      num_eliminated_supernodes_(0),
      num_aggressive_absorptions_(0) {
  QUOTIENT_START_TIMER(timers_, kSetup);

  graph_data_.signed_supernode_sizes.Resize(num_vertices_, 1);
  graph_data_.dense_supernode.size = 0;
  graph_data_.dense_supernode.principal_member = -1;

  // Initialize the counts for the adjacency lists.
  graph_data_.element_offsets.Resize(num_vertices_);
  graph_data_.element_sizes.Resize(num_vertices_, 0);
  graph_data_.adjacency_list_sizes.Resize(num_vertices_, 0);
  for (const MatrixEntry<Field>& entry : entries) {
    if (entry.row != entry.column) {
      ++graph_data_.adjacency_list_sizes[entry.row];
    }
  }

  Int num_edges = ConvertEdgeCountsIntoOffsets();

  // Pack the edges. We leave extra room for packing element structures.
  // TODO(Jack Poulson): Make this coefficient configurable.
  const float kNumEdgesRatio = 0.2f;
  const Int extra_element_space = kNumEdgesRatio * num_edges + num_vertices_;
  graph_data_.lists.Resize(num_edges + extra_element_space);
  num_edges = 0;
  for (const MatrixEntry<Field>& entry : entries) {
    const Int source = entry.row;
    if (graph_data_.signed_supernode_sizes[source] && entry.column != source) {
      graph_data_.lists[num_edges++] = entry.column;
    }
  }
  graph_data_.offset = num_edges;

  InitializeDegreeAndHashLists();
  InitializeNodeFlags();

  QUOTIENT_STOP_TIMER(timers_, kSetup);
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
  pivot_ = degrees_and_hashes_.lists.FindMinimalIndex(
      control_.force_minimal_pivot_indices);
  return pivot_;
}

inline Int QuotientGraph::NumAggressiveAbsorptions() const QUOTIENT_NOEXCEPT {
  return num_aggressive_absorptions_;
}

inline Int QuotientGraph::NumDense() const QUOTIENT_NOEXCEPT {
  return graph_data_.dense_supernode.size;
}

inline Int QuotientGraph::NumVertices() const QUOTIENT_NOEXCEPT {
  return num_vertices_;
}

inline Int QuotientGraph::NumEliminatedVertices() const QUOTIENT_NOEXCEPT {
  return num_eliminated_vertices_;
}

inline Int QuotientGraph::NumHashBucketCollisions() const QUOTIENT_NOEXCEPT {
  return degrees_and_hashes_.num_bucket_collisions;
}

inline Int QuotientGraph::NumHashCollisions() const QUOTIENT_NOEXCEPT {
  return degrees_and_hashes_.num_collisions;
}

inline Int QuotientGraph::SupernodeSize(Int i) const QUOTIENT_NOEXCEPT {
  if (graph_data_.signed_supernode_sizes[i]) {
    // This is a traditional supernode (but it might be eliminated).
    return std::abs(graph_data_.signed_supernode_sizes[i]);
  } else if (graph_data_.ActiveSupernode(i)) {
    // This is a dense node.
    return 1;
  } else {
    // This is a merged node.
    return 0;
  }
}

inline const Buffer<Int> QuotientGraph::Element(Int i) const QUOTIENT_NOEXCEPT {
  const Int element_size = graph_data_.ElementSize(i);
  const Int* element_data = graph_data_.ElementData(i);

  Buffer<Int> element_copy(element_size);
  for (Int j = 0; j < element_size; ++j) {
    element_copy[j] = element_data[j];
  }

  return element_copy;
}

inline Buffer<Int> QuotientGraph::ElementList(Int i) const QUOTIENT_NOEXCEPT {
  const Int num_elements = graph_data_.ElementListSize(i);
  const Int* element_list = graph_data_.ElementList(i);

  Buffer<Int> element_list_copy(num_elements);
  for (Int k = 0; k < num_elements; ++k) {
    element_list_copy[k] = element_list[k];
  }

  return element_list_copy;
}

inline void QuotientGraph::ComputePivotStructure() QUOTIENT_NOEXCEPT {
  QUOTIENT_START_TIMER(timers_, kComputePivotStructure);
  const Int num_elements = graph_data_.ElementListSize(pivot_);
  const Int num_adjacencies = graph_data_.adjacency_list_sizes[pivot_];
  const Int pivot_supernode_size = graph_data_.signed_supernode_sizes[pivot_];

  // Allocate space for the element using an upper-bound on the size
  // (note that, because of supernodes, this is *not* the degree).
  Int element_size_bound = num_adjacencies;
  for (Int k = 0; k < num_elements; ++k) {
    const Int* element_list = graph_data_.ElementList(pivot_);
    const Int element = element_list[k];
    const Int element_size = graph_data_.ElementSize(element);
    QUOTIENT_ASSERT(element_size > 0, "Non-positive element size of " +
                                          std::to_string(element_size));
    element_size_bound += element_size - 1;
  }
  if (graph_data_.offset + element_size_bound > Int(graph_data_.lists.Size())) {
#ifdef QUOTIENT_DEBUG
    std::cout << "Repacking with offset: " << graph_data_.offset
              << ", element_size_bound: " << element_size_bound
              << ", lists.Size(): " << graph_data_.lists.Size() << std::endl;
#endif
    graph_data_.Pack();
  }
  const Int* element_list = graph_data_.ElementList(pivot_);
  const Int* adjacency_list = graph_data_.AdjacencyList(pivot_);

  graph_data_.element_offsets[pivot_] = graph_data_.offset;
  Int* pivot_data = graph_data_.ElementData(pivot_);

  // Negate the signed supernode size of the pivot.
  graph_data_.signed_supernode_sizes[pivot_] = -pivot_supernode_size;

  // Push the supervariables in the pivot adjacency list into the structure.
  Int pivot_size = 0;
  Int pivot_degree = 0;
  for (Int k = 0; k < num_adjacencies; ++k) {
    const Int i = adjacency_list[k];
    const Int supernode_size = graph_data_.signed_supernode_sizes[i];
    QUOTIENT_ASSERT(supernode_size >= 0,
                    "An element was in the adjacency list.");
    if (!supernode_size) {
      continue;
    }
    QUOTIENT_ASSERT(graph_data_.ActiveSupernode(i),
                    "An absorbed element was in the adjacency list.");

    pivot_degree += supernode_size;
    graph_data_.signed_supernode_sizes[i] = -supernode_size;

    // Store index 'i' into the pivot element.
    pivot_data[pivot_size++] = i;
  }
  graph_data_.adjacency_list_sizes[pivot_] = 0;

  // Push the unique supervariables in the patterns of the pivot element list
  // into the current structure.
  for (Int k = 0; k < num_elements; ++k) {
    const Int element = element_list[k];
    QUOTIENT_ASSERT(graph_data_.ActiveSupernode(element),
                    "Used an absorbed element in pivot structure.");

    const Int element_size = graph_data_.ElementSize(element);
    const Int* element_data = graph_data_.ElementData(element);
    for (Int j = 0; j < element_size; ++j) {
      const Int index = element_data[j];
      const Int supernode_size = graph_data_.signed_supernode_sizes[index];
      // While no eliminated supernodes should appear in unabsorbed elements,
      // we have (temporarily) flipped the signs of the supernode sizes of
      // the members we have already added to this pivot's structure.
      if (supernode_size <= 0) {
        continue;
      }
      QUOTIENT_ASSERT(graph_data_.ActiveSupernode(index),
                      "An absorbed element was in an element.");

      pivot_degree += supernode_size;
      graph_data_.signed_supernode_sizes[index] = -supernode_size;

      // Store index 'index' into the pivot element.
      pivot_data[pivot_size++] = index;
    }
#ifdef QUOTIENT_DEBUG
    Int degree = 0;
    for (Int j = 0; j < element_size; ++j) {
      const Int index = element_data[j];
      const Int supernode_size = -graph_data_.signed_supernode_sizes[index];
      QUOTIENT_ASSERT(supernode_size >= 0,
                      "Flipped supernode size was expected to be positive.");
      degree += supernode_size;
    }
    QUOTIENT_ASSERT(degree == degrees_and_hashes_.lists.degrees[element],
                    "Degree did not match its cached value.");
#endif

    // Absorb this element into the pivot.
    degrees_and_hashes_.lists.degrees[element] -= pivot_supernode_size;
    graph_data_.SetParent(element, pivot_);
    graph_data_.ElementSize(element) = 0;
    node_flags_.flags[element] = 0;
  }

  // Set the size of the pivot element.
  graph_data_.ElementSize(pivot_) = pivot_size;
  QUOTIENT_ASSERT(
      graph_data_.offset + pivot_size <= Int(graph_data_.lists.Size()),
      "Packed beyond end of element indices.");

  degrees_and_hashes_.lists.degrees[pivot_] = pivot_degree;
  node_flags_.max_degree = std::max(node_flags_.max_degree, pivot_degree);

  QUOTIENT_STOP_TIMER(timers_, kComputePivotStructure);
}

inline Int QuotientGraph::NumPivotElements() const QUOTIENT_NOEXCEPT {
  return graph_data_.ElementListSize(pivot_);
}

inline Int QuotientGraph::NumPivotDegreeUpdates() const QUOTIENT_NOEXCEPT {
  return graph_data_.ElementSize(pivot_);
}

inline Int QuotientGraph::NumPivotDegreeUpdatesWithMultipleElements() const
    QUOTIENT_NOEXCEPT {
  Int num_multi_updates = 0;
  const Int pivot_size = graph_data_.ElementSize(pivot_);
  const Int* pivot_data = graph_data_.ElementData(pivot_);
  for (Int j = 0; j < pivot_size; ++j) {
    const Int i = pivot_data[j];
    if (graph_data_.ElementListSize(i) > 2) {
      ++num_multi_updates;
    }
  }
  return num_multi_updates;
}

inline Int QuotientGraph::NumPivotCholeskyNonzeros() const QUOTIENT_NOEXCEPT {
  const Int pivot_size = std::abs(graph_data_.signed_supernode_sizes[pivot_]);
  const Int structure_size = degrees_and_hashes_.lists.degrees[pivot_] +
                             graph_data_.dense_supernode.size;
  const Int diag_block_nonzeros = (pivot_size * (pivot_size + 1)) / 2;
  const Int subdiagonal_nonzeros = structure_size * pivot_size;
  return diag_block_nonzeros + subdiagonal_nonzeros;
}

inline double QuotientGraph::NumPivotCholeskyFlops() const QUOTIENT_NOEXCEPT {
  const Int pivot_size = std::abs(graph_data_.signed_supernode_sizes[pivot_]);
  const Int structure_size = degrees_and_hashes_.lists.degrees[pivot_] +
                             graph_data_.dense_supernode.size;
  const double diag_block_flops = std::pow(1. * pivot_size, 3.) / 3.;
  const double schur_complement_flops =
      std::pow(1. * structure_size, 2.) * pivot_size;
  return diag_block_flops + schur_complement_flops;
}

inline void QuotientGraph::PackCountAndHashAdjacencies(
    Int i, Int num_elements, Int* degree, UInt* hash) QUOTIENT_NOEXCEPT {
  Int degree_new = *degree;
  UInt hash_new = *hash;

  const Int old_num_elements = graph_data_.ElementListSize(i);
  QUOTIENT_ASSERT(num_elements <= old_num_elements,
                  "Packing adjacencies after where they begin.");
  const Int num_adjacencies = graph_data_.adjacency_list_sizes[i];
  const Int offset = graph_data_.element_offsets[i];
  const Int* old_adjacency_list = &graph_data_.lists[offset + old_num_elements];
  Int* adjacency_list = &graph_data_.lists[offset + num_elements];

  Int num_packed = 0;
  for (Int k = 0; k < num_adjacencies; ++k) {
    const Int j = old_adjacency_list[k];
    // We must filter out: non-principal variables, eliminated supernodes,
    // the pivot supernode, and all members of the pivot structure.
    //
    // Due to our (temporary) negation of the supernode sizes of the pivot and
    // members of the pivot structure, all four cases can be handled by
    // demanding a positive (signed) supernode size.
    const Int supernode_size = graph_data_.signed_supernode_sizes[j];
    if (supernode_size > 0) {
      degree_new += supernode_size;
      QUOTIENT_HASH_COMBINE(hash_new, j);
      adjacency_list[num_packed++] = j;
    }
  }
  graph_data_.ElementListSize(i) = num_elements;
  graph_data_.adjacency_list_sizes[i] = num_packed;
  *degree = degree_new;
  *hash = hash_new;
}

inline void QuotientGraph::InsertPivotElement(Int i) QUOTIENT_NOEXCEPT {
  Int& num_elements = graph_data_.ElementListSize(i);
  const Int element_list_beg = graph_data_.element_offsets[i];
  const Int element_list_end = element_list_beg + num_elements;
  const Int adjacency_list_end =
      element_list_end + graph_data_.adjacency_list_sizes[i];
#ifdef QUOTIENT_DEBUG
  if (i == num_vertices_ - 1) {
    if (adjacency_list_end == Int(graph_data_.lists.Size())) {
      std::cerr << "Adjacency list ran off the end of the array." << std::endl;
    }
  } else {
    if (adjacency_list_end == graph_data_.element_offsets[i + 1]) {
      std::cerr << "Adjacency list overlapped with next element list."
                << std::endl;
    }
  }
#endif
  if (graph_data_.adjacency_list_sizes[i]) {
    // Move the first adjacency to the back.
    graph_data_.lists[adjacency_list_end] = graph_data_.lists[element_list_end];
  }
  graph_data_.lists[element_list_end] = pivot_;
  ++num_elements;
}

inline std::pair<Int, UInt> QuotientGraph::ExactEmptyDegreeAndHash(Int i)
    QUOTIENT_NOEXCEPT {
  const Int num_elements = graph_data_.ElementListSize(i);
  Int degree = 0;
  UInt hash = 0;

  // We should only have one member of the element list, 'pivot'.
  QUOTIENT_ASSERT(
      num_elements == 0,
      "The element list was falsely assumed to have a single entry.");

  // Add |L_p \ supernode(i)|.
  const Int supernode_size = -graph_data_.signed_supernode_sizes[i];
  QUOTIENT_ASSERT(supernode_size > 0,
                  "The flipped supernode size should have been positive.");
  degree += degrees_and_hashes_.lists.degrees[pivot_] - supernode_size;
  QUOTIENT_HASH_COMBINE(hash, pivot_);

  // Add the cardinality of A_i \ supernode(i), where 'i' is the principal
  // variable.
  PackCountAndHashAdjacencies(i, num_elements, &degree, &hash);
  InsertPivotElement(i);

  return std::make_pair(degree, hash);
}

inline std::pair<Int, UInt> QuotientGraph::ExactSingleDegreeAndHash(Int i)
    QUOTIENT_NOEXCEPT {
  Int num_elements = graph_data_.ElementListSize(i);
  const Int* element_list = graph_data_.ElementList(i);
  Int degree = 0;
  UInt hash = 0;

  // There should be exactly two members in the element list, and the second
  // should be the pivot.
  QUOTIENT_ASSERT(num_elements == 1,
                  "The element list should have had exactly one members.");
  const Int element = element_list[0];
  QUOTIENT_ASSERT(element != pivot_,
                  "The element list should have contained a non-pivot.");

  // Add |L_p \ supernode(i)|.
  const Int supernode_size = -graph_data_.signed_supernode_sizes[i];
  QUOTIENT_ASSERT(supernode_size > 0,
                  "The flipped supernode size should have been positive.");
  degree += degrees_and_hashes_.lists.degrees[pivot_] - supernode_size;
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
  Int num_elements = graph_data_.ElementListSize(i);
  Int* element_list = graph_data_.ElementList(i);
  Int degree = 0;
  UInt hash = 0;
  const Int shift = node_flags_.shift;

  // Add on the number of unique entries in the structures of the element lists
  // that are outside supernode(i).
  Int num_packed = 0;
  for (Int k = 0; k < num_elements; ++k) {
    const Int element = element_list[k];
    if (!node_flags_.flags[element]) {
      continue;
    }

    element_list[num_packed++] = element;
    QUOTIENT_HASH_COMBINE(hash, element);

    const Int element_size = graph_data_.ElementSize(element);
    const Int* element_data = graph_data_.ElementData(element);
    for (Int index = 0; index < element_size; ++index) {
      const Int j = element_data[index];

      // Unabsorbed elements should not have any eliminated members of
      // in their element. Thus, we can take the absolute value of the signed
      // supernode size (since some members might be in the pivot structure).
      if (node_flags_.flags[j] == shift || i == j) {
        continue;
      }
      degree += std::abs(graph_data_.signed_supernode_sizes[j]);
      node_flags_.flags[j] = shift;
    }
  }
  num_elements = num_packed;

  // Handle the pivot (but don't pack it yet).
  QUOTIENT_HASH_COMBINE(hash, pivot_);
  const Int pivot_size = graph_data_.ElementSize(pivot_);
  const Int* pivot_data = graph_data_.ElementData(pivot_);
  for (Int index = 0; index < pivot_size; ++index) {
    const Int j = pivot_data[index];

    // Unabsorbed elements should not have any eliminated members of
    // in their element. Thus, we can take the absolute value of the signed
    // supernode size (since some members might be in the pivot structure).
    if (node_flags_.flags[j] == shift || i == j) {
      continue;
    }
    degree += std::abs(graph_data_.signed_supernode_sizes[j]);
    node_flags_.flags[j] = shift;
  }

  // Clear the mask.
  for (Int k = 0; k < num_elements; ++k) {
    const Int element = element_list[k];
    QUOTIENT_ASSERT(node_flags_.flags[element],
                    "Flag was unset after packing.");
    const Int element_size = graph_data_.ElementSize(element);
    const Int* element_data = graph_data_.ElementData(element);
    for (Int index = 0; index < element_size; ++index) {
      const Int j = element_data[index];
      node_flags_.flags[j] = shift - 1;
    }
  }
  for (Int index = 0; index < pivot_size; ++index) {
    const Int j = pivot_data[index];
    node_flags_.flags[j] = shift - 1;
  }

  // Add the cardinality of A_i \ supernode(i), where 'i' is the principal
  // variable.
  PackCountAndHashAdjacencies(i, num_elements, &degree, &hash);
  InsertPivotElement(i);

  return std::make_pair(degree, hash);
}

inline void QuotientGraph::ExactDegreesAndHashes() QUOTIENT_NOEXCEPT {
  const Int pivot_size = graph_data_.ElementSize(pivot_);
  const Int* pivot_data = graph_data_.ElementData(pivot_);
  for (Int index = 0; index < pivot_size; ++index) {
    const Int i = pivot_data[index];

    std::pair<Int, UInt> degree_and_hash;
    const Int num_elements = graph_data_.ElementListSize(i);
    if (num_elements == 0) {
      degree_and_hash = ExactEmptyDegreeAndHash(i);
    } else if (num_elements == 1) {
      degree_and_hash = ExactSingleDegreeAndHash(i);
    } else {
      degree_and_hash = ExactGenericDegreeAndHash(i);
    }
    QUOTIENT_ASSERT(
        degree_and_hash.first >= 0,
        "Tried to set a degree of " + std::to_string(degree_and_hash.first));
    degrees_and_hashes_.lists.degrees[i] = degree_and_hash.first;
    degrees_and_hashes_.lists.SetHash(i, degree_and_hash.second);
  }
}

inline void QuotientGraph::AmestoyDegreesAndHashes() QUOTIENT_NOEXCEPT {
  const Int shift = node_flags_.shift;
  const Int pivot_degree = degrees_and_hashes_.lists.degrees[pivot_];
  const Int pivot_size = graph_data_.ElementSize(pivot_);
  const Int* pivot_data = graph_data_.ElementData(pivot_);
  const Int num_vertices_left = num_vertices_ - num_eliminated_vertices_;

  for (Int index = 0; index < pivot_size; ++index) {
    const Int i = pivot_data[index];
    Int degree = 0;
    UInt hash = 0;

    const Int supernode_size = -graph_data_.signed_supernode_sizes[i];
    QUOTIENT_ASSERT(supernode_size > 0,
                    "The negated supernode size was expected to be positive.");
    const Int bound0 = num_vertices_left - supernode_size;

    // Note that this usage of 'external' refers to |L_p \ supernode(i)| and not
    // |L_e \ L_p|, as is the case for 'external_degrees'.
    const Int external_pivot_degree = pivot_degree - supernode_size;
    QUOTIENT_ASSERT(external_pivot_degree >= 0,
                    "Encountered a negative external pivot degree");

    const Int old_degree = degrees_and_hashes_.lists.degrees[i];
    const Int bound1 = old_degree + external_pivot_degree;

    // bound_2 = |A_i \ supernode(i)| + |L_p \ supernode(i)| +
    //           \sum_{e in E_i \ {p}} |L_e \ L_p|.
    degree += external_pivot_degree;
    Int num_packed = 0;
    Int num_elements = graph_data_.ElementListSize(i);
    Int* element_list = graph_data_.ElementList(i);
    for (Int k = 0; k < num_elements; ++k) {
      const Int element = element_list[k];
      QUOTIENT_ASSERT(element != pivot_, "Iterated over pivot element.");
      const Int flag = node_flags_.flags[element];
      if (flag) {
        element_list[num_packed++] = element;

        const Int external_degree = flag - shift;
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

    QUOTIENT_ASSERT(degree >= 0,
                    "Tried to set a degree of " + std::to_string(degree));
    degrees_and_hashes_.lists.degrees[i] = degree;
    degrees_and_hashes_.lists.SetHash(i, hash);
  }
}

inline void QuotientGraph::AshcraftDegreesAndHashes() QUOTIENT_NOEXCEPT {
  const Int pivot_size = graph_data_.ElementSize(pivot_);
  const Int* pivot_data = graph_data_.ElementData(pivot_);
  for (Int index = 0; index < pivot_size; ++index) {
    const Int i = pivot_data[index];

    std::pair<Int, UInt> degree_and_hash;
    // Note that the list size is one *before* adding the pivot.
    if (graph_data_.ElementListSize(i) == 1) {
      degree_and_hash = ExactSingleDegreeAndHash(i);
    } else {
      degree_and_hash = GilbertDegreeAndHash(i);
    }
    QUOTIENT_ASSERT(
        degree_and_hash.first >= 0,
        "Tried to set a degree of " + std::to_string(degree_and_hash.first));
    degrees_and_hashes_.lists.degrees[i] = degree_and_hash.first;
    degrees_and_hashes_.lists.SetHash(i, degree_and_hash.second);
  }
}

inline std::pair<Int, UInt> QuotientGraph::GilbertDegreeAndHash(Int i)
    QUOTIENT_NOEXCEPT {
  Int degree = 0;
  UInt hash = 0;

  const Int supernode_size = -graph_data_.signed_supernode_sizes[i];
  QUOTIENT_ASSERT(supernode_size > 0,
                  "The negated supernode size was expected to be positive.");

  Int num_packed = 0;
  Int num_elements = graph_data_.ElementListSize(i);
  Int* element_list = graph_data_.ElementList(i);
  for (Int k = 0; k < num_elements; ++k) {
    const Int element = element_list[k];
    if (!node_flags_.flags[element]) {
      continue;
    }
    element_list[num_packed++] = element;

    QUOTIENT_ASSERT(graph_data_.ActiveSupernode(element),
                    "Used absorbed element in Gilbert degree update.");
    QUOTIENT_ASSERT(
        degrees_and_hashes_.lists.degrees[element] - supernode_size >= 0,
        "Negative Gilbert degree update.");
    degree += degrees_and_hashes_.lists.degrees[element] - supernode_size;
    QUOTIENT_HASH_COMBINE(hash, element);
  }
  num_elements = num_packed;
  degree += degrees_and_hashes_.lists.degrees[pivot_] - supernode_size;
  QUOTIENT_HASH_COMBINE(hash, pivot_);

  PackCountAndHashAdjacencies(i, num_elements, &degree, &hash);
  InsertPivotElement(i);

  const Int num_vertices_left = num_vertices_ - num_eliminated_vertices_;
  const Int bound = num_vertices_left - supernode_size;
  degree = std::min(degree, bound);

  return std::make_pair(degree, hash);
}

inline void QuotientGraph::GilbertDegreesAndHashes() QUOTIENT_NOEXCEPT {
  const Int num_vertices_left = num_vertices_ - num_eliminated_vertices_;
  const Int pivot_degree = degrees_and_hashes_.lists.degrees[pivot_];
  const Int pivot_size = graph_data_.ElementSize(pivot_);
  const Int* pivot_data = graph_data_.ElementData(pivot_);
  for (Int index = 0; index < pivot_size; ++index) {
    const Int i = pivot_data[index];
    Int degree = 0;
    UInt hash = 0;

    const Int supernode_size = -graph_data_.signed_supernode_sizes[i];
    QUOTIENT_ASSERT(supernode_size > 0,
                    "The negated supernode size was expected to be positive.");

    Int num_packed = 0;
    Int num_elements = graph_data_.ElementListSize(i);
    Int* element_list = graph_data_.ElementList(i);
    for (Int k = 0; k < num_elements; ++k) {
      const Int element = element_list[k];
      if (!node_flags_.flags[element]) {
        continue;
      }
      element_list[num_packed++] = element;

      QUOTIENT_ASSERT(graph_data_.ActiveSupernode(element),
                      "Used absorbed element in Gilbert degree update.");
      QUOTIENT_ASSERT(
          degrees_and_hashes_.lists.degrees[element] - supernode_size >= 0,
          "Negative Gilbert degree update.");
      degree += degrees_and_hashes_.lists.degrees[element] - supernode_size;
      QUOTIENT_HASH_COMBINE(hash, element);
    }
    num_elements = num_packed;
    degree += pivot_degree - supernode_size;
    QUOTIENT_HASH_COMBINE(hash, pivot_);

    PackCountAndHashAdjacencies(i, num_elements, &degree, &hash);
    InsertPivotElement(i);

    const Int bound = num_vertices_left - supernode_size;
    degree = std::min(degree, bound);

    QUOTIENT_ASSERT(degree >= 0,
                    "Tried to set a degree of " + std::to_string(degree));
    degrees_and_hashes_.lists.degrees[i] = degree;
    degrees_and_hashes_.lists.SetHash(i, hash);
  }
}

inline void QuotientGraph::ComputeDegreesAndHashes() QUOTIENT_NOEXCEPT {
  QUOTIENT_START_TIMER(timers_, kComputeDegreesAndHashes);

  // Remove the old degree values from the linked lists.
  // (We return many of them when finalizing the pivot.)
  const Int pivot_size = graph_data_.ElementSize(pivot_);
  const Int* pivot_data = graph_data_.ElementData(pivot_);
  for (Int index = 0; index < pivot_size; ++index) {
    const Int i = pivot_data[index];
    degrees_and_hashes_.lists.RemoveDegree(i);
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
  QUOTIENT_STOP_TIMER(timers_, kComputeDegreesAndHashes);
}

inline bool QuotientGraph::StructuralVariablesAreQuotientIndistinguishable(
    Int i, Int j) const QUOTIENT_NOEXCEPT {
  QUOTIENT_ASSERT(
      i != j,
      "Explicitly tested for equivalence of a supernode with itself...");
  const Int num_elements_i = graph_data_.ElementListSize(i);
  const Int* element_list_i = graph_data_.ElementList(i);
  const Int* element_list_j = graph_data_.ElementList(j);

  const Int num_adjacencies_j = graph_data_.adjacency_list_sizes[j];
  const Int* adjacency_list_j = graph_data_.AdjacencyList(j);

#ifdef QUOTIENT_DEBUG
  for (Int k = 0; k < num_elements_i; ++k) {
    const Int element = element_list_i[k];
    QUOTIENT_ASSERT(
        graph_data_.ActiveSupernode(element),
        "Absorbed element was in element list during absorption test.");
  }
  const Int num_elements_j = graph_data_.ElementListSize(j);
  for (Int k = 0; k < num_elements_j; ++k) {
    const Int element = element_list_j[k];
    QUOTIENT_ASSERT(
        graph_data_.ActiveSupernode(element),
        "Absorbed element was in element list during absorption test.");
  }
  const Int num_adjacencies_i = graph_data_.adjacency_list_sizes[i];
  const Int* adjacency_list_i = graph_data_.AdjacencyList(i);
  for (Int k = 0; k < num_adjacencies_i; ++k) {
    const Int index = adjacency_list_i[k];
    QUOTIENT_ASSERT(
        graph_data_.signed_supernode_sizes[index] > 0,
        "Non-positive supernode size in adj list during absorption test.");
  }
  for (Int k = 0; k < num_adjacencies_j; ++k) {
    const Int index = adjacency_list_j[k];
    QUOTIENT_ASSERT(
        graph_data_.signed_supernode_sizes[index] > 0,
        "Non-positive supernode size in adj list during absorption test.");
  }
#endif

  // It is assumed that the number of elements and adjacencies is equivalent.
  QUOTIENT_ASSERT(num_elements_i == num_elements_j,
                  "Mismatched number of elements.");
  QUOTIENT_ASSERT(num_adjacencies_i == num_adjacencies_j,
                  "Mismatched number of adjacencies.");

  // Check if A_i = A_j by ensuring that all members of A_j are in A_i (by
  // checking that the indices are currently flagged).
  const Int shift = node_flags_.shift;
  for (Int k = 0; k < num_adjacencies_j; ++k) {
    const Int index = adjacency_list_j[k];
    if (node_flags_.flags[index] != shift) {
      return false;
    }
  }

  // Check if E_i = E_j.
  for (Int k = 0; k < num_elements_i; ++k) {
    if (element_list_i[k] != element_list_j[k]) {
      // TODO(Jack Poulson): Assert that the k'th member of the element_list_i
      // is not contained anywhere in element_list_j.
      return false;
    }
  }

  return true;
}

inline void QuotientGraph::MergeVariables() QUOTIENT_NOEXCEPT {
  QUOTIENT_START_TIMER(timers_, kMergeVariables);
  const Int pivot_size = graph_data_.ElementSize(pivot_);
  const Int* pivot_data = graph_data_.ElementData(pivot_);
  const Int shift = node_flags_.shift;

  // Add the hashes into the hash lists.
  for (Int i_index = 0; i_index < pivot_size; ++i_index) {
    const Int i = pivot_data[i_index];
    const UInt hash = degrees_and_hashes_.lists.Hash(i);
    const Int bucket = hash % num_vertices_;
    degrees_and_hashes_.lists.AddHash(i, hash, bucket);
  }

  Int num_merges = 0;
  const Buffer<Int>& next_member = degrees_and_hashes_.lists.next_member;
  for (Int i_index = 0; i_index < pivot_size; ++i_index) {
    Int i = pivot_data[i_index];
    const UInt hash = degrees_and_hashes_.lists.Hash(i);
    const Int bucket = hash % num_vertices_;
    const Int head = degrees_and_hashes_.lists.HashBucketHead(bucket);
    // We are not the head of the bucket.
    if (i != head) {
      continue;
    }

    // Test the unique pairs in the bucket.
    for (; next_member[i] != -1; i = next_member[i]) {
      if (!graph_data_.signed_supernode_sizes[i]) {
        continue;
      }
      QUOTIENT_ASSERT(graph_data_.signed_supernode_sizes[i] < 0,
                      "Supernode size should have been temporarily negative.");
      const UInt i_hash = degrees_and_hashes_.lists.Hash(i);
      bool scattered_adjacencies = false;

      const Int num_adjacencies_i = graph_data_.adjacency_list_sizes[i];
      const Int* adjacency_list_i = graph_data_.AdjacencyList(i);

      for (Int j = next_member[i]; j != -1; j = next_member[j]) {
        if (!graph_data_.signed_supernode_sizes[j]) {
          continue;
        }
        QUOTIENT_ASSERT(
            graph_data_.signed_supernode_sizes[j] < 0,
            "Supernode size should have been temporarily negative.");
        const UInt j_hash = degrees_and_hashes_.lists.Hash(j);
        if (i_hash != j_hash ||
            graph_data_.ElementListSize(i) != graph_data_.ElementListSize(j) ||
            num_adjacencies_i != graph_data_.adjacency_list_sizes[j]) {
          ++degrees_and_hashes_.num_bucket_collisions;
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
          for (Int k = 0; k < num_adjacencies_i; ++k) {
            const Int index = adjacency_list_i[k];
            node_flags_.flags[index] = shift;
          }
          scattered_adjacencies = true;
        }

        if (StructuralVariablesAreQuotientIndistinguishable(i, j)) {
          ++num_merges;
          const Int absorbed_size = -graph_data_.signed_supernode_sizes[j];
          QUOTIENT_ASSERT(absorbed_size > 0,
                          "Absorbed size should have been positive.");

          // Merge [i] -> [j] (and i becomes the principal member).
          //
          // Recall that the supernode sizes of principal variables in the
          // pivot structure have been negated.
          degrees_and_hashes_.lists.degrees[i] -= absorbed_size;
          graph_data_.signed_supernode_sizes[i] -= absorbed_size;
          graph_data_.SetParent(j, i);
          graph_data_.ElementListSize(j) = 0;
          graph_data_.adjacency_list_sizes[j] = 0;
          graph_data_.signed_supernode_sizes[j] = 0;
        } else {
          ++degrees_and_hashes_.num_collisions;
          ++degrees_and_hashes_.num_bucket_collisions;
        }
      }
      if (scattered_adjacencies) {
        // Reset the flagged adjacencies to the shift minus one.
        for (Int k = 0; k < num_adjacencies_i; ++k) {
          const Int index = adjacency_list_i[k];
          node_flags_.flags[index] = shift - 1;
        }
      }
    }
    degrees_and_hashes_.lists.ClearHashBucket(bucket);
  }

#ifdef QUOTIENT_DEBUG
  for (Int i_index = 0; i_index < pivot_size; ++i_index) {
    Int i = pivot_data[i_index];
    const UInt hash = degrees_and_hashes_.lists.Hash(i);
    const Int bucket = hash % num_vertices_;
    const Int index = degrees_and_hashes_.lists.HashBucketHead(bucket);
    if (index != -1) {
      std::cerr << "Did not clear head for bucket " << bucket << std::endl;
    }
  }
#endif

  QUOTIENT_STOP_TIMER(timers_, kMergeVariables);
}

inline void QuotientGraph::FinalizePivot() QUOTIENT_NOEXCEPT {
  QUOTIENT_START_TIMER(timers_, kFinalizePivot);
  ResetExternalDegrees();

  const Int supernode_size = -graph_data_.signed_supernode_sizes[pivot_];
  QUOTIENT_ASSERT(supernode_size > 0,
                  "The supernode size was assumed positive.");

  // Unflip the signs of supernode sizes for the members of the pivot structure
  // and simultaneously remove the non-principal members.
  Int num_packed = 0;
  const Int pivot_size = graph_data_.ElementSize(pivot_);
  Int* pivot_data = graph_data_.ElementData(pivot_);
  for (Int index = 0; index < pivot_size; ++index) {
    const Int i = pivot_data[index];
    QUOTIENT_ASSERT(
        graph_data_.signed_supernode_sizes[i] <= 0,
        "A member of pivot element had a positive signed supernode size.");
    if (graph_data_.signed_supernode_sizes[i] < 0) {
      const Int degree = degrees_and_hashes_.lists.degrees[i];
      QUOTIENT_ASSERT(degree >= 0,
                      "Negative degree of " + std::to_string(degree) +
                          " at pivot structure index " + std::to_string(index) +
                          " (" + std::to_string(i) + ") for pivot " +
                          std::to_string(pivot_));
      degrees_and_hashes_.lists.AddDegree(i, degree);
      graph_data_.signed_supernode_sizes[i] =
          -graph_data_.signed_supernode_sizes[i];
      pivot_data[num_packed++] = i;
    }
  }
  graph_data_.ElementSize(pivot_) = num_packed;
  graph_data_.offset += num_packed;

  ++num_eliminated_supernodes_;
  num_eliminated_vertices_ += supernode_size;

  QUOTIENT_STOP_TIMER(timers_, kFinalizePivot);
}

inline void QuotientGraph::ComputePostorder(Buffer<Int>* postorder)
    QUOTIENT_NOEXCEPT {
  QUOTIENT_START_TIMER(timers_, kComputePostorder);

  // Fill the supernode non-principal member lists.
  // Also, ensure that the root nodes are explicitly marked as such.
  Buffer<Int> nonprincipal_offsets(num_vertices_ + 1, 0);
  Int num_nonprincipal_members = 0;
  for (Int i = 0; i < num_vertices_; ++i) {
    nonprincipal_offsets[i] = num_nonprincipal_members;
    const Int supernode_size = -graph_data_.signed_supernode_sizes[i];
    QUOTIENT_ASSERT(supernode_size >= 0, "Supernode size was negative.");
    if (supernode_size) {
      num_nonprincipal_members += supernode_size - 1;
      if (graph_data_.element_offsets[i] >= 0) {
        graph_data_.SetParent(i, -1);
      }
    }
  }
  nonprincipal_offsets[num_vertices_] = num_nonprincipal_members;
  auto offsets_copy = nonprincipal_offsets;
  Buffer<Int> nonprincipal_members(num_nonprincipal_members);
  for (Int i = 0; i < num_vertices_; ++i) {
    const Int supernode_size = -graph_data_.signed_supernode_sizes[i];
    if (!supernode_size) {
      Int principal = i;
      while (!graph_data_.signed_supernode_sizes[principal]) {
        QUOTIENT_ASSERT(
            !graph_data_.ActiveSupernode(principal),
            "Active supernode while computing supernode principal.");
        principal = graph_data_.Parent(principal);
      }
      nonprincipal_members[offsets_copy[principal]++] = i;
    }
  }

  // Reconstruct the child links from the parent links in a contiguous array
  // (similar to the CSR format) by first counting the number of children of
  // each node.

  Buffer<Int> child_offsets;
  {
    Buffer<Int> num_children(num_vertices_, 0);
    for (Int i = 0; i < num_vertices_; ++i) {
      if (!graph_data_.signed_supernode_sizes[i]) {
        continue;
      }
      const Int parent = graph_data_.Parent(i);
      if (parent >= 0) {
        ++num_children[parent];
      }
    }
    OffsetScan(num_children, &child_offsets);
  }
  Buffer<Int> children(num_eliminated_supernodes_);
  {
    auto offsets_copy = child_offsets;
    for (Int i = 0; i < num_vertices_; ++i) {
      if (!graph_data_.signed_supernode_sizes[i]) {
        continue;
      }
      const Int parent = graph_data_.Parent(i);
      if (parent >= 0) {
        children[offsets_copy[parent]++] = i;
      }
    }
  }

  // Scan for the roots and launch a pe-order traversal on each of them.
  postorder->Resize(num_vertices_);
  Int* iter = postorder->begin();
  for (Int i = 0; i < num_vertices_; ++i) {
    if (!graph_data_.ActiveSupernode(i)) {
      // This element was absorbed into another element.
      continue;
    }
    iter = PreorderTree(i, nonprincipal_members, nonprincipal_offsets, children,
                        child_offsets, iter);
  }

  // Reverse the preordering (to form a postordering) in-place.
  std::reverse(postorder->begin(), iter);

  if (graph_data_.dense_supernode.size) {
    // Pack the dense columns into the back.
    for (Int i = 0; i < num_vertices_; ++i) {
      if (!graph_data_.signed_supernode_sizes[i] &&
          graph_data_.ActiveSupernode(i)) {
        *iter = i;
        ++iter;
      }
    }
  }

  QUOTIENT_ASSERT(iter == postorder->end(),
                  "Postorder had incorrect final offset.");
  QUOTIENT_STOP_TIMER(timers_, kComputePostorder);
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
  permuted_supernode_sizes->Clear();
  permuted_supernode_sizes->Resize(num_eliminated_supernodes_);

  Int i_perm = 0;
  for (Int index = 0; index < num_eliminated_supernodes_; ++index) {
    const Int i = inverse_permutation[i_perm];
    QUOTIENT_ASSERT(graph_data_.signed_supernode_sizes[i] < 0,
                    "Supernode size was negative.");
    const Int supernode_size = -graph_data_.signed_supernode_sizes[i];
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
    if (graph_data_.signed_supernode_sizes[i]) {
      // 'i' is a principal variable of an absorbed element.
      QUOTIENT_ASSERT(graph_data_.signed_supernode_sizes[i] < 0,
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
    Buffer<Int>* permuted_graph_data_parents) const QUOTIENT_NOEXCEPT {
  permuted_graph_data_parents->Clear();
  permuted_graph_data_parents->Resize(num_eliminated_supernodes_);
  for (Int index = 0; index < num_vertices_; ++index) {
    if (!graph_data_.signed_supernode_sizes[index]) {
      continue;
    }
    const Int parent = graph_data_.Parent(index);

    const Int permuted_index = permutation[index];
    const Int permuted_parent = parent >= 0 ? permutation[parent] : -1;

    const Int permuted_supernode = permuted_member_to_supernode[permuted_index];
    const Int permuted_parent_supernode =
        permuted_parent >= 0 ? permuted_member_to_supernode[permuted_parent]
                             : -1;

    (*permuted_graph_data_parents)[permuted_supernode] =
        permuted_parent_supernode;
  }
}

inline void QuotientGraph::ExternalDegrees() QUOTIENT_NOEXCEPT {
  QUOTIENT_START_TIMER(timers_, kExternalDegrees);
  const Int shift = node_flags_.shift;
  const bool aggressive_absorption = control_.aggressive_absorption;

  const Int pivot_size = graph_data_.ElementSize(pivot_);
  const Int* pivot_data = graph_data_.ElementData(pivot_);
  for (Int index = 0; index < pivot_size; ++index) {
    const Int i = pivot_data[index];
    const Int supernode_size = -graph_data_.signed_supernode_sizes[i];
    QUOTIENT_ASSERT(supernode_size > 0,
                    "supernode " + std::to_string(i) +
                        " had non-positive signed size "
                        "when computing external element sizes");
    const Int shift_minus_supernode_size = shift - supernode_size;

    const Int num_elements = graph_data_.ElementListSize(i);
    const Int* element_list = graph_data_.ElementList(i);
    for (Int k = 0; k < num_elements; ++k) {
      const Int element = element_list[k];
      Int& shifted_external_degree = node_flags_.flags[element];

      if (shifted_external_degree >= shift) {
        QUOTIENT_ASSERT(graph_data_.ActiveSupernode(element),
                        "Tried to subtract using an absorbed element.");
        QUOTIENT_ASSERT(
            shifted_external_degree != shift,
            "Shifted external size was equal to shift before subtracting");
        shifted_external_degree -= supernode_size;
      } else if (shifted_external_degree) {
        QUOTIENT_ASSERT(graph_data_.ActiveSupernode(element),
                        "Tried to subtract using an absorbed element.");
        shifted_external_degree = degrees_and_hashes_.lists.degrees[element] +
                                  shift_minus_supernode_size;
      }
      QUOTIENT_ASSERT(
          !shifted_external_degree || shifted_external_degree >= shift,
          "Computed negative external element degree.");
      if (aggressive_absorption && shifted_external_degree == shift) {
        ++num_aggressive_absorptions_;
        shifted_external_degree = 0;
        graph_data_.SetParent(element, pivot_);
        graph_data_.ElementSize(element) = 0;
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
  for (std::size_t i = 0; i < node_flags_.flags.Size(); ++i) {
    if (node_flags_.flags[i]) {
      node_flags_.flags[i] = 1;
    }
  }
}

inline void QuotientGraph::CombineDenseNodes() QUOTIENT_NOEXCEPT {
  if (!graph_data_.dense_supernode.size) {
    return;
  }

  // Absorb all of the dense nodes into the principal member.
  for (Int i = 0; i < num_vertices_; ++i) {
    if (!graph_data_.signed_supernode_sizes[i] &&
        graph_data_.ActiveSupernode(i)) {
      if (i == graph_data_.dense_supernode.principal_member) {
        graph_data_.signed_supernode_sizes[i] =
            -graph_data_.dense_supernode.size;
      } else {
        graph_data_.SetParent(i, graph_data_.dense_supernode.principal_member);
        graph_data_.signed_supernode_sizes[i] = 0;
      }
    }
  }

  // Point the non-dense elements currently marked as roots to the dense
  // supernode.
  for (Int i = 0; i < num_vertices_; ++i) {
    if (graph_data_.signed_supernode_sizes[i] &&
        i != graph_data_.dense_supernode.principal_member &&
        graph_data_.ActiveSupernode(i)) {
      graph_data_.SetParent(i, graph_data_.dense_supernode.principal_member);
    }
  }
}

inline const MinimumDegreeControl& QuotientGraph::Control() const
    QUOTIENT_NOEXCEPT {
  return control_;
}

inline std::vector<std::pair<std::string, double>>
QuotientGraph::ComponentSeconds() const QUOTIENT_NOEXCEPT {
  std::vector<std::pair<std::string, double>> times;
#ifdef QUOTIENT_ENABLE_TIMERS
  for (const std::pair<std::string, Timer>& pairing : timers_) {
    times.emplace_back(pairing.first, pairing.second.TotalSeconds());
  }
#endif
  return times;
}

}  // namespace quotient

#endif  // ifndef QUOTIENT_QUOTIENT_GRAPH_IMPL_H_
