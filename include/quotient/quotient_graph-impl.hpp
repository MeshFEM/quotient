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

#include "quotient/config.hpp"
#include "quotient/coordinate_graph.hpp"
#include "quotient/macros.hpp"
#include "quotient/quotient_graph.hpp"

namespace quotient {

#ifdef QUOTIENT_ENABLE_TIMERS
#  define QUOTIENT_START_TIMER(timer, name) timer[name].Start()
#  define QUOTIENT_STOP_TIMER(timer, name) timer[name].Stop()
#else
#  define QUOTIENT_START_TIMER(timer, name)
#  define QUOTIENT_STOP_TIMER(timer, name)
#endif

#ifdef QUOTIENT_ENABLE_TIMERS
static constexpr char kSetup[] = "Setup";
static constexpr char kComputePivotStructure[] = "ComputePivotStructure";
static constexpr char kExternalDegrees[] = "ExternalDegrees";
static constexpr char kComputeDegrees[] = "ComputeDegrees";
static constexpr char kMergeVariables[] = "MergeVariables";
static constexpr char kFinalizePivot[] = "FinalizePivot";
#endif

template<typename T>
void PrintVector(const std::vector<T>& vec, const std::string& msg) {
  std::cout << msg << ": ";
  for (UInt i = 0; i < vec.size(); ++i) {
    std::cout << vec[i] << " ";
  }
  std::cout << "\n";
}


inline QuotientGraph::QuotientGraph(
    const CoordinateGraph& graph,
    const MinimumDegreeControl& control)
: num_original_vertices_(graph.NumSources()),
  num_eliminated_vertices_(0),
  control_(control),
  max_degree_(0),
  external_degree_shift_(0),
  max_shift_value_(std::numeric_limits<Int>::max() - graph.NumSources()),
  num_hash_bucket_collisions_(0),
  num_hash_collisions_(0),
  num_aggressive_absorptions_(0),
  num_dense_(0) {
  QUOTIENT_START_TIMER(timers_, kSetup);

  // Initialize the supernodes as simple.
  signed_supernode_sizes_.resize(num_original_vertices_, 1);
  elimination_order_.reserve(num_original_vertices_);
  parents_.resize(num_original_vertices_, -1);
  next_index_.resize(num_original_vertices_, -1);
  tail_index_.resize(num_original_vertices_);
  std::iota(tail_index_.begin(), tail_index_.end(), 0);

  // Initialize the counts for the adjacency lists.
  element_list_offsets_.resize(num_original_vertices_, 0);
  element_list_sizes_.resize(num_original_vertices_, 0);
  adjacency_list_sizes_.resize(num_original_vertices_, 0);
  const std::vector<GraphEdge>& edges = graph.Edges();
  for (Int source = 0; source < num_original_vertices_; ++source) {
    const Int source_edge_offset = graph.SourceEdgeOffset(source);
    const Int next_source_edge_offset = graph.SourceEdgeOffset(source + 1);
    for (Int k = source_edge_offset; k < next_source_edge_offset; ++k) {
      const GraphEdge& edge = edges[k];
      if (edge.second != source) {
        ++adjacency_list_sizes_[source];
      }
    }
  }

  // Every row with at least this many non-diagonal nonzeros will be treated
  // as dense and moved to the back of the postordering.
  const Int dense_threshold =
      std::max(1.f * control_.min_dense_threshold,
          control_.dense_sqrt_multiple *
              std::sqrt(1.f * num_original_vertices_));

  // Convert the counts into an offset scan.
  Int num_edges = 0;
  for (Int source = 0; source < num_original_vertices_; ++source) {
    element_list_offsets_[source] = num_edges;
    if (adjacency_list_sizes_[source] >= dense_threshold) {
      ++num_dense_;
      ++num_eliminated_vertices_;
      // We will denote this source as dense by setting its supernode size to
      // 0 and keeping its parent as -1. A merged variable also has a supernode
      // size of 0, but its parent is a valid vertex index.
      signed_supernode_sizes_[source] = 0;
      adjacency_list_sizes_[source] = 0;
    }
    num_edges += adjacency_list_sizes_[source];
  }
  element_list_offsets_[num_original_vertices_] = num_edges;
#ifdef QUOTIENT_DEBUG
  if (num_dense_) {
    std::cout << "Eliminated " << num_dense_ << " dense rows." << std::endl;
  }
#endif

  // Pack the edges.
  num_edges = 0;
  element_and_adjacency_lists_.reserve(num_edges);
  for (Int source = 0; source < num_original_vertices_; ++source) {
    if (!signed_supernode_sizes_[source]) {
      // Skip the dense row.
      continue;
    }
    const Int source_edge_offset = graph.SourceEdgeOffset(source);
    const Int next_source_edge_offset = graph.SourceEdgeOffset(source + 1);
    for (Int k = source_edge_offset; k < next_source_edge_offset; ++k) {
      const GraphEdge& edge = edges[k];
      if (edge.second != source) {
        element_and_adjacency_lists_.push_back(edge.second);
      }
    }
  }

  // Initialize the degree lists.
  degree_lists_.degrees.resize(num_original_vertices_, 0);
  degree_lists_.heads.resize(num_original_vertices_ - 1, -1);
  degree_lists_.next_member.resize(num_original_vertices_, -1);
  degree_lists_.last_member.resize(num_original_vertices_, -1);
  for (Int source = 0; source < num_original_vertices_; ++source) {
    if (!signed_supernode_sizes_[source]) {
      // Skip the dense row.
      continue;
    }
    const Int degree = adjacency_list_sizes_[source];
    degree_lists_.AddDegree(source, degree);
  }

  // Initialize the hash lists.
  hash_lists_.buckets.resize(num_original_vertices_);
  hash_lists_.hashes.resize(num_original_vertices_);
  hash_lists_.heads.resize(num_original_vertices_, -1);
  hash_lists_.next_member.resize(num_original_vertices_, -1);

  // Trivially initialize the lower-triangular nonzero structures.
  if (control_.store_structures) {
    structures_.resize(num_original_vertices_);
  }
  elements_.resize(num_original_vertices_);

  // The absorbed elements will be maintained as 0, the unabsorbed will be
  // initialized as 1, and the shift will be initialized as 2.
  external_degree_shift_ = 2;
  node_flags_.resize(num_original_vertices_, 1);

  QUOTIENT_STOP_TIMER(timers_, kSetup);
}

inline const std::vector<Int>& QuotientGraph::EliminationOrder() const {
  return elimination_order_;
}

inline Int QuotientGraph::FindAndProcessPivot() {
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

inline Int QuotientGraph::GetNextPivot() {
  pivot_ = degree_lists_.FindMinimalIndex(control_.force_minimal_pivot_indices);
  return pivot_;
}

inline Int QuotientGraph::NumAggressiveAbsorptions() const {
  return num_aggressive_absorptions_;
}

inline Int QuotientGraph::NumDense() const {
  return num_dense_;
}

inline Int QuotientGraph::NumOriginalVertices() const {
  return num_original_vertices_;
}

inline Int QuotientGraph::NumEliminatedVertices() const {
  return num_eliminated_vertices_;
}

inline Int QuotientGraph::NumHashBucketCollisions() const {
  return num_hash_bucket_collisions_;
}

inline Int QuotientGraph::NumHashCollisions() const {
  return num_hash_collisions_;
}

inline void QuotientGraph::Print() const {
  for (Int i = 0; i < num_original_vertices_; ++i) {
    if (!signed_supernode_sizes_[i]) {
      continue;
    }
    std::cout << "Supernode " << i << "\n";
    const std::vector<Int> supernode = FormSupernode(i);
    PrintVector(supernode, "  members");
    PrintVector(elements_[i], "  elements");

    std::cout << "  element list: ";
    const Int element_list_beg = element_list_offsets_[i];
    const Int element_list_end = element_list_beg + element_list_sizes_[i];
    for (Int k = element_list_beg; k < element_list_end; ++k) {
      std::cout << element_and_adjacency_lists_[k] << " ";
    }
    std::cout << "\n";

    std::cout << "  adjacency list: ";
    const Int adjacency_list_end = element_list_end + adjacency_list_sizes_[i];
    for (Int k = element_list_end; k < adjacency_list_end; ++k) {
      std::cout << element_and_adjacency_lists_[k] << " ";
    }
    std::cout << "\n";
    std::cout << "\n";
  }
}

inline std::vector<Int> QuotientGraph::FormSupernode(Int i) const {
  const Int supernode_size = SupernodeSize(i);
  std::vector<Int> supernode;
  AppendSupernode(i, supernode_size, &supernode);
  return supernode;
}

inline Int QuotientGraph::SupernodeSize(Int i) const {
  if (signed_supernode_sizes_[i]) {
    // This is a traditional supernode (but it might be eliminated).
    return std::abs(signed_supernode_sizes_[i]);
  } else if (parents_[i] == -1) {
    // This is a dense node.
    return 1;
  } else {
    // This is a merged node.
    return 0;
  }
}

inline const std::vector<Int>& QuotientGraph::Element(Int i) const {
  return elements_[i];
}

inline std::vector<Int> QuotientGraph::ElementList(Int i) const {
  const Int num_elements = element_list_sizes_[i];
  const Int element_list_beg = element_list_offsets_[i];
  const Int element_list_end = element_list_beg + num_elements;
  std::vector<Int> element_list;
  element_list.reserve(num_elements);
  for (Int k = element_list_beg; k < element_list_end; ++k) {
    element_list.push_back(element_and_adjacency_lists_[k]);
  }
  return element_list;
}

inline void QuotientGraph::FormEliminatedStructures(
    std::vector<std::vector<Int>>* eliminated_structures) const {
  eliminated_structures->resize(elimination_order_.size());
  for (UInt index = 0; index < elimination_order_.size(); ++index) {
    std::vector<Int>& eliminated_structure = (*eliminated_structures)[index];
    eliminated_structure = structures_[elimination_order_[index]];
  }
}

inline void QuotientGraph::ComputePivotStructure() {
  QUOTIENT_START_TIMER(timers_, kComputePivotStructure);
  QUOTIENT_ASSERT(elements_[pivot_].empty(), "Chose a pivot more than once.");
  const Int element_list_beg = element_list_offsets_[pivot_];
  const Int element_list_end = element_list_beg + element_list_sizes_[pivot_];
  const Int adjacency_list_end =
      element_list_end + adjacency_list_sizes_[pivot_];
  const Int pivot_supernode_size = signed_supernode_sizes_[pivot_];
  std::vector<Int>& pivot_element = elements_[pivot_];

  // Allocate space for the element using an upper-bound on the size
  // (note that, because of supernodes, this is *not* the degree).
  Int element_size_bound = adjacency_list_sizes_[pivot_];
  for (Int k = element_list_beg; k < element_list_end; ++k) {
    const Int element = element_and_adjacency_lists_[k];
    element_size_bound += elements_[element].size() - 1;
  }
  pivot_element.reserve(element_size_bound);

  // Negate the signed supernode size of the pivot.
  signed_supernode_sizes_[pivot_] = -pivot_supernode_size;

  // Push the supervariables in the pivot adjacency list into the structure.
  Int pivot_degree = 0;
  for (Int k = element_list_end; k < adjacency_list_end; ++k) {
    const Int i = element_and_adjacency_lists_[k];
    const Int supernode_size = signed_supernode_sizes_[i];
    QUOTIENT_ASSERT(supernode_size >= 0,
        "An element was in the adjacency list.");
    if (!supernode_size) {
      continue;
    }
    QUOTIENT_ASSERT(parents_[i] == -1,
        "An absorbed element was in the adjacency list.");

    pivot_degree += supernode_size;
    signed_supernode_sizes_[i] = -supernode_size;
    pivot_element.push_back(i);
  }
  adjacency_list_sizes_[pivot_] = 0;

  // Push the unique supervariables in the patterns of the pivot element list
  // into the current structure.
  for (Int k = element_list_beg; k < element_list_end; ++k) {
    const Int element = element_and_adjacency_lists_[k];
    QUOTIENT_ASSERT(parents_[element] == -1,
        "Used an absorbed element in pivot structure.");

    for (const Int& index : elements_[element]) {
      const Int supernode_size = signed_supernode_sizes_[index];
      // While no eliminated supernodes should appear in unabsorbed elements,
      // we have (temporarily) flipped the signs of the supernode sizes of
      // the members we have already added to this pivot's structure.
      if (supernode_size <= 0) {
        continue;
      }
      QUOTIENT_ASSERT(parents_[index] == -1,
          "An absorbed element was in an element.");

      pivot_degree += supernode_size;
      signed_supernode_sizes_[index] = -supernode_size;
      pivot_element.push_back(index);
    }
#ifdef QUOTIENT_DEBUG
    Int degree = 0;
    for (const Int& index : elements_[element]) {
      const Int supernode_size = -signed_supernode_sizes_[index];
      QUOTIENT_ASSERT(supernode_size >= 0,
          "Flipped supernode size was expected to be positive.");
      degree += supernode_size;
    }
    QUOTIENT_ASSERT(degree == degree_lists_.degrees[element],
        "Degree did not match its cached value.");
#endif

    // Absorb this element into the pivot.
    parents_[element] = pivot_;
    node_flags_[element] = 0;
    SwapClearVector(&elements_[element]);
    element_list_sizes_[element] = 0;
  }

  degree_lists_.degrees[pivot_] = pivot_degree;
  max_degree_ = std::max(max_degree_, pivot_degree);

  if (control_.store_structures) {
    for (const Int& i : elements_[pivot_]) {
      const Int supernode_size = -signed_supernode_sizes_[i];
      AppendSupernode(i, supernode_size, &structures_[pivot_]);
    }
  }

  QUOTIENT_STOP_TIMER(timers_, kComputePivotStructure);
}

inline Int QuotientGraph::NumPivotElements() const {
  return element_list_sizes_[pivot_];
}

inline Int QuotientGraph::NumPivotDegreeUpdates() const {
  return elements_[pivot_].size();
}

inline Int QuotientGraph::NumPivotDegreeUpdatesWithMultipleElements() const {
  Int num_multi_updates = 0;
  for (const Int& i : elements_[pivot_]) {
    if (element_list_sizes_[i] > 2) {
      ++num_multi_updates;
    }
  }
  return num_multi_updates;
}

inline Int QuotientGraph::NumPivotCholeskyNonzeros() const {
  const Int pivot_size = std::abs(signed_supernode_sizes_[pivot_]);
  const Int structure_size = degree_lists_.degrees[pivot_] + num_dense_;
  const Int diag_block_nonzeros = (pivot_size * (pivot_size + 1)) / 2;
  const Int subdiagonal_nonzeros = structure_size * pivot_size;
  return diag_block_nonzeros + subdiagonal_nonzeros;
}

inline double QuotientGraph::NumPivotCholeskyFlops() const {
  const Int pivot_size = std::abs(signed_supernode_sizes_[pivot_]);
  const Int structure_size = degree_lists_.degrees[pivot_] + num_dense_;
  const double diag_block_flops =  std::pow(1. * pivot_size, 3.) / 3.;
  const double schur_complement_flops =
      std::pow(1. * structure_size, 2.) * pivot_size;
  return diag_block_flops + schur_complement_flops;
}

inline void QuotientGraph::PackCountAndHashAdjacencies(
    Int i, Int num_elements, Int* degree, UInt* hash) {
  Int degree_new = *degree;
  UInt hash_new = *hash;

  Int num_packed = 0;
  const Int pack_index = element_list_offsets_[i] + num_elements;
  const Int adjacency_list_beg =
      element_list_offsets_[i] + element_list_sizes_[i];
  QUOTIENT_ASSERT(pack_index <= adjacency_list_beg,
      "Packing adjacencies after where they begin.");
  const Int adjacency_list_end = adjacency_list_beg + adjacency_list_sizes_[i];
  for (Int k = adjacency_list_beg; k < adjacency_list_end; ++k) {
    const Int j = element_and_adjacency_lists_[k];
    // We must filter out: non-principal variables, eliminated supernodes,
    // the pivot supernode, and all members of the pivot structure.
    //
    // Due to our (temporary) negation of the supernode sizes of the pivot and
    // members of the pivot structure, all four cases can be handled by
    // demanding a positive (signed) supernode size.
    const Int supernode_size = signed_supernode_sizes_[j];
    if (supernode_size > 0) {
      degree_new += supernode_size;
      QUOTIENT_HASH_COMBINE(hash_new, j);
      element_and_adjacency_lists_[pack_index + num_packed++] = j;
    }
  }
  element_list_sizes_[i] = num_elements;
  adjacency_list_sizes_[i] = num_packed;
  *degree = degree_new;
  *hash = hash_new;
}

inline void QuotientGraph::InsertPivotElement(Int i) {
  Int& num_elements = element_list_sizes_[i];
  const Int element_list_beg = element_list_offsets_[i];
  const Int element_list_end = element_list_beg + num_elements;
  const Int adjacency_list_end = element_list_end + adjacency_list_sizes_[i];
#ifdef QUOTIENT_DEBUG
  if (i == num_original_vertices_ - 1) {
    if (adjacency_list_end == Int(element_and_adjacency_lists_.size())) {
      std::cerr << "Adjacency list ran off the end of the array." << std::endl;
    }
  } else {
    if (adjacency_list_end == element_list_offsets_[i + 1]) {
      std::cerr << "Adjacency list overlapped with next element list."
                << std::endl;
    }
  }
#endif
  if (adjacency_list_sizes_[i]) {
    // Move the first adjacency to the back.
    element_and_adjacency_lists_[adjacency_list_end] =
        element_and_adjacency_lists_[element_list_end];
  }
  element_and_adjacency_lists_[element_list_end] = pivot_;
  ++num_elements;
}

inline std::pair<Int, UInt>
QuotientGraph::ExactEmptyDegreeAndHash(Int i) {
  const Int num_elements = element_list_sizes_[i];
  Int degree = 0;
  UInt hash = 0;

  // We should only have one member of the element list, 'pivot'.
  QUOTIENT_ASSERT(num_elements == 0,
      "The element list was falsely assumed to have a single entry.");

  // Add |L_p \ supernode(i)|.
  const Int supernode_size = -signed_supernode_sizes_[i];
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

inline std::pair<Int, UInt>
QuotientGraph::ExactSingleDegreeAndHash(Int i) {
  const Int offset = element_list_offsets_[i];
  Int num_elements = element_list_sizes_[i];
  Int degree = 0;
  UInt hash = 0;

  // There should be exactly two members in the element list, and the second
  // should be the pivot.
  QUOTIENT_ASSERT(num_elements == 1,
      "The element list should have had exactly one members.");
  const Int& element = element_and_adjacency_lists_[offset];
  QUOTIENT_ASSERT(element != pivot_,
      "The element list should have contained a non-pivot.");

  // Add |L_p \ supernode(i)|.
  const Int supernode_size = -signed_supernode_sizes_[i];
  QUOTIENT_ASSERT(supernode_size > 0,
      "The flipped supernode size should have been positive.");
  degree += degree_lists_.degrees[pivot_] - supernode_size;
  QUOTIENT_HASH_COMBINE(hash, pivot_);

  Int& shifted_external_degree = node_flags_[element];
  if (!shifted_external_degree) {
    --num_elements;
  } else {
    // Add |L_e \ L_p|.
    degree += shifted_external_degree - external_degree_shift_;
    QUOTIENT_HASH_COMBINE(hash, element);
  }

  // Add the cardinality of A_i \ supernode(i), where 'i' is the principal
  // variable.
  PackCountAndHashAdjacencies(i, num_elements, &degree, &hash);
  InsertPivotElement(i);

  return std::make_pair(degree, hash);
}

inline std::pair<Int, UInt>
QuotientGraph::ExactGenericDegreeAndHash(Int i) {
  const Int offset = element_list_offsets_[i];
  Int num_elements = element_list_sizes_[i];
  Int degree = 0;
  UInt hash = 0;
  const Int shift = external_degree_shift_;

  // Add on the number of unique entries in the structures of the element lists
  // that are outside supernode(i).
  Int num_packed = 0;
  for (Int k = 0; k < num_elements; ++k) {
    const Int element = element_and_adjacency_lists_[offset + k];
    if (!node_flags_[element]) {
      continue;
    }

    element_and_adjacency_lists_[offset + num_packed++] = element;
    QUOTIENT_HASH_COMBINE(hash, element);

    for (const Int& j : elements_[element]) {
      // Unabsorbed elements should not have any eliminated members of
      // in their element. Thus, we can take the absolute value of the signed
      // supernode size (since some members might be in the pivot structure).
      if (node_flags_[j] == shift || i == j) {
        continue;
      }
      degree += std::abs(signed_supernode_sizes_[j]);
      node_flags_[j] = shift;
    }
  }
  num_elements = num_packed;

  // Handle the pivot (but don't pack it yet).
  QUOTIENT_HASH_COMBINE(hash, pivot_);
  for (const Int& j : elements_[pivot_]) {
    // Unabsorbed elements should not have any eliminated members of
    // in their element. Thus, we can take the absolute value of the signed
    // supernode size (since some members might be in the pivot structure).
    if (node_flags_[j] == shift || i == j) {
      continue;
    }
    degree += std::abs(signed_supernode_sizes_[j]);
    node_flags_[j] = shift;
  }

  // Clear the mask.
  for (Int k = 0; k < num_elements; ++k) {
    const Int element = element_and_adjacency_lists_[offset + k];
    if (!node_flags_[element]) {
      continue;
    }
    for (const Int& j : elements_[element]) {
      node_flags_[j] = shift - 1;
    }
  }
  for (const Int& j : elements_[pivot_]) {
    node_flags_[j] = shift - 1;
  }

  // Add the cardinality of A_i \ supernode(i), where 'i' is the principal
  // variable.
  PackCountAndHashAdjacencies(i, num_elements, &degree, &hash);
  InsertPivotElement(i);

  return std::make_pair(degree, hash);
}

inline void QuotientGraph::ExactDegreesAndHashes() {
  for (const Int& i : elements_[pivot_]) {
    std::pair<Int, UInt> degree_and_hash;
    const Int num_elements = element_list_sizes_[i];
    if (num_elements == 0) {
      degree_and_hash = ExactEmptyDegreeAndHash(i);
    } else if (num_elements == 1) {
      degree_and_hash = ExactSingleDegreeAndHash(i);
    } else {
      degree_and_hash = ExactGenericDegreeAndHash(i);
    }
    degree_lists_.degrees[i] = degree_and_hash.first;
    hash_lists_.hashes[i] = degree_and_hash.second;
  }
}

inline void QuotientGraph::AmestoyDegreesAndHashes() {
  const Int shift = external_degree_shift_;
  const Int pivot_degree = degree_lists_.degrees[pivot_];
  const Int num_vertices_left =
      num_original_vertices_ - num_eliminated_vertices_;
  for (const Int& i : elements_[pivot_]) {
    Int degree = 0;
    UInt hash = 0;

    const Int supernode_size = -signed_supernode_sizes_[i];
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
    Int num_elements = element_list_sizes_[i];
    const Int offset = element_list_offsets_[i];
    for (Int k = offset; k < offset + num_elements; ++k) {
      const Int element = element_and_adjacency_lists_[k];
      QUOTIENT_ASSERT(element != pivot_, "Iterated over pivot element.");
      if (node_flags_[element]){
        QUOTIENT_ASSERT(node_flags_[element],
            "Ran into an absorbed element in the external degree calculation.");
        element_and_adjacency_lists_[offset + num_packed++] = element;

        const Int external_degree = node_flags_[element] - shift;
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
    hash_lists_.hashes[i] = hash;
  }
}

inline void QuotientGraph::AshcraftDegreesAndHashes() {
  for (const Int& i : elements_[pivot_]) {
    std::pair<Int, UInt> degree_and_hash;
    // Note that the list size is one *before* adding the pivot.
    if (element_list_sizes_[i] == 1) {
      degree_and_hash = ExactSingleDegreeAndHash(i);
    } else {
      degree_and_hash = GilbertDegreeAndHash(i);
    }
    degree_lists_.degrees[i] = degree_and_hash.first;
    hash_lists_.hashes[i] = degree_and_hash.second;
  }
}

inline std::pair<Int, UInt> QuotientGraph::GilbertDegreeAndHash(Int i) {
  Int degree = 0;
  UInt hash = 0;

  const Int supernode_size = -signed_supernode_sizes_[i];
  QUOTIENT_ASSERT(supernode_size > 0,
      "The negated supernode size was expected to be positive.");

  Int num_packed = 0;
  const Int offset = element_list_offsets_[i];
  Int num_elements = element_list_sizes_[i];
  for (Int k = 0; k < num_elements; ++k) {
    const Int element = element_and_adjacency_lists_[offset + k];
    if (!node_flags_[element]) {
      continue;
    }
    element_and_adjacency_lists_[offset + num_packed++] = element;

    QUOTIENT_ASSERT(parents_[element] == -1,
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

  const Int num_vertices_left =
      num_original_vertices_ - num_eliminated_vertices_;
  const Int bound = num_vertices_left - supernode_size;
  degree = std::min(degree, bound);

  return std::make_pair(degree, hash);
}

inline void QuotientGraph::GilbertDegreesAndHashes() {
  const Int num_vertices_left =
      num_original_vertices_ - num_eliminated_vertices_;
  const Int pivot_degree = degree_lists_.degrees[pivot_];
  for (const Int& i : elements_[pivot_]) {
    Int degree = 0;
    UInt hash = 0;

    const Int supernode_size = -signed_supernode_sizes_[i];
    QUOTIENT_ASSERT(supernode_size > 0,
        "The negated supernode size was expected to be positive.");

    Int num_packed = 0;
    Int num_elements = element_list_sizes_[i];
    const Int offset = element_list_offsets_[i];
    for (Int k = 0; k < num_elements; ++k) {
      const Int element = element_and_adjacency_lists_[offset + k];
      if (!node_flags_[element]) {
        continue;
      }
      element_and_adjacency_lists_[offset + num_packed++] = element;

      QUOTIENT_ASSERT(parents_[element] == -1,
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
    hash_lists_.hashes[i] = hash;
  }
}

inline void QuotientGraph::ComputeDegreesAndHashes() {
  QUOTIENT_START_TIMER(timers_, kComputeDegrees);

  // Remove the old values from the linked lists.
  for (const Int& i : elements_[pivot_]) {
    degree_lists_.RemoveDegree(i);
  }

  switch(control_.degree_type) {
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
    Int i, Int j) const {
  QUOTIENT_ASSERT(i != j,
      "Explicitly tested for equivalence of a supernode with itself...");
  const Int offset_i = element_list_offsets_[i];
  const Int offset_j = element_list_offsets_[j];
  const Int num_elements_i = element_list_sizes_[i];
  const Int num_elements_j = element_list_sizes_[j];
  const Int num_adjacencies_i = adjacency_list_sizes_[i];
  const Int num_adjacencies_j = adjacency_list_sizes_[j];

#ifdef QUOTIENT_DEBUG
  for (Int k = 0; k < num_elements_i; ++k) {
    const Int element = element_and_adjacency_lists_[offset_i + k];
    QUOTIENT_ASSERT(parents_[element] == -1,
        "Absorbed element was in element list during absorption test.");
  }
  for (Int k = 0; k < num_elements_j; ++k) {
    const Int element = element_and_adjacency_lists_[offset_j + k];
    QUOTIENT_ASSERT(parents_[element] == -1,
        "Absorbed element was in element list during absorption test.");
  }
  for (Int k = 0; k < num_adjacencies_i; ++k) {
    const Int index =
        element_and_adjacency_lists_[offset_i + num_elements_i + k];
    QUOTIENT_ASSERT(signed_supernode_sizes_[index] > 0,
        "Non-positive supernode size in adj list during absorption test.");
  }
  for (Int k = 0; k < num_adjacencies_j; ++k) {
    const Int index =
        element_and_adjacency_lists_[offset_j + num_elements_j + k];
    QUOTIENT_ASSERT(signed_supernode_sizes_[index] > 0,
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
    if (element_and_adjacency_lists_[offset_i + k] !=
        element_and_adjacency_lists_[offset_j + k]) {
      return false;
    }
  }

  // Check if A_i = A_j by ensuring that all members of A_j are in A_i (by
  // checking that the indices are currently flagged).
  const Int shift = external_degree_shift_;
  const Int adjacency_list_j_beg = offset_j + num_elements_j;
  const Int adjacency_list_j_end = adjacency_list_j_beg + num_adjacencies_j;
  for (Int k = adjacency_list_j_beg; k < adjacency_list_j_end; ++k) {
    const Int index = element_and_adjacency_lists_[k];
    if (node_flags_[index] != shift) {
      return false;
    }
  }

  return true;
}

inline void QuotientGraph::MergeVariables() {
  QUOTIENT_START_TIMER(timers_, kMergeVariables);
  const std::vector<Int>& pivot_element = elements_[pivot_];
  const Int supernodal_struct_size = pivot_element.size();
  const Int shift = external_degree_shift_;

  // Add the hashes into the hash lists.
  for (Int i_index = 0; i_index < supernodal_struct_size; ++i_index) {
    const Int i = pivot_element[i_index];
    const UInt hash = hash_lists_.hashes[i];
    const Int bucket = hash % num_original_vertices_;
    hash_lists_.AddHash(i, hash, bucket);
  }

  Int num_merges = 0;
  const std::vector<Int>& next_member = hash_lists_.next_member;
  for (Int i_index = 0; i_index < supernodal_struct_size; ++i_index) {
    Int i = pivot_element[i_index];
    const Int bucket = hash_lists_.buckets[i];
    const Int head = hash_lists_.heads[bucket];
    // We are not the head of the bucket.
    if (i != head) {
      continue;
    }

    // Test the unique pairs in the bucket.
    for (; next_member[i] != -1; i = next_member[i]) {
      if (!signed_supernode_sizes_[i]) {
        continue;
      }
      QUOTIENT_ASSERT(signed_supernode_sizes_[i] < 0,
          "Supernode size should have been temporarily negative.");
      const UInt i_hash = hash_lists_.hashes[i];
      bool scattered_adjacencies = false;
      const Int adjacency_list_i_beg =
          element_list_offsets_[i] + element_list_sizes_[i];
      const Int adjacency_list_i_end =
          adjacency_list_i_beg + adjacency_list_sizes_[i];
      for (Int j = next_member[i]; j != -1; j = next_member[j]) {
        if (!signed_supernode_sizes_[j]) { 
          continue;
        }
        QUOTIENT_ASSERT(signed_supernode_sizes_[j] < 0,
            "Supernode size should have been temporarily negative.");
        const UInt j_hash = hash_lists_.hashes[j];
        if (i_hash != j_hash) {
          ++num_hash_bucket_collisions_;
          continue;
        }
        if (!scattered_adjacencies) {
          // Temporarily flag the adjacencies of supervariable i by setting
          // 'node_flags_' to the current shift value in said indices. This
          // is representative because only elements should hold this value.
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
          for (Int k = adjacency_list_i_beg; k < adjacency_list_i_end; ++k) {
            const Int index = element_and_adjacency_lists_[k];
            node_flags_[index] = shift;
          }
          scattered_adjacencies = true;
        }

        if (StructuralVariablesAreQuotientIndistinguishable(i, j)) {
          ++num_merges;
          const Int absorbed_size = -signed_supernode_sizes_[j];
          QUOTIENT_ASSERT(absorbed_size > 0,
              "Absorbed size should have been positive.");
 
          //#pragma omp critical
          const Int old_degree = degree_lists_.degrees[i];
          const Int degree = old_degree - absorbed_size;
          degree_lists_.UpdateDegree(i, degree);
          degree_lists_.RemoveDegree(j);

          // Merge [i] -> [j] (and i becomes the principal member).
          parents_[j] = i;
          next_index_[tail_index_[i]] = j;
          tail_index_[i] = tail_index_[j];
          signed_supernode_sizes_[i] -= absorbed_size;  // Recall it is negated.
          signed_supernode_sizes_[j] = 0;

          element_list_sizes_[j] = 0;
          adjacency_list_sizes_[j] = 0;
        } else {
          ++num_hash_bucket_collisions_;
          ++num_hash_collisions_;
        }
      }
      if (scattered_adjacencies) {
        // Reset the flagged adjacencies to the shift minus one.
        for (Int k = adjacency_list_i_beg; k < adjacency_list_i_end; ++k) {
          const Int index = element_and_adjacency_lists_[k];
          node_flags_[index] = shift - 1;
        }
      }
    }
    hash_lists_.ClearBucket(bucket);
  }

#ifdef QUOTIENT_DEBUG
  for (Int index = 0; index < num_original_vertices_; ++index) {
    if (hash_lists_.heads[index] != -1) {
      std::cerr << "Did not clear head for bucket " << index << std::endl;
    }
  }
#endif

  QUOTIENT_STOP_TIMER(timers_, kMergeVariables);
}

inline void QuotientGraph::FinalizePivot() {
  QUOTIENT_START_TIMER(timers_, kFinalizePivot);
  ResetExternalDegrees();

  const Int supernode_size = -signed_supernode_sizes_[pivot_];
  QUOTIENT_ASSERT(supernode_size > 0,
      "The supernode size was assumed positive.");

  // Since this supervariable is being eliminated, it needs to be implicitly
  // removed from all elements containing it. This is accomplished through
  // decrementing their element sizes and marking the supernode as eliminated
  // by making its sign negative.
  const Int offset = element_list_offsets_[pivot_];
  const Int num_elements = element_list_sizes_[pivot_];
  for (Int k = offset; k < offset + num_elements; ++k) {
    const Int element = element_and_adjacency_lists_[k];
    degree_lists_.degrees[element] -= supernode_size;
  }

  // Unflip the signs of supernode sizes for the members of the pivot structure
  // and simultaneously remove the non-principal members.
  Int num_packed = 0;
  for (const Int& i : elements_[pivot_]) {
    QUOTIENT_ASSERT(signed_supernode_sizes_[i] <= 0,
        "A member of pivot element had a positive signed supernode size.");
    if (signed_supernode_sizes_[i] < 0) {
      signed_supernode_sizes_[i] = -signed_supernode_sizes_[i];
      elements_[pivot_][num_packed++] = i;
    }
  }
  elements_[pivot_].resize(num_packed);

  elimination_order_.push_back(pivot_);
  num_eliminated_vertices_ += supernode_size;

  QUOTIENT_STOP_TIMER(timers_, kFinalizePivot);
}

inline void QuotientGraph::AppendSupernode(
    Int principal_member, Int supernode_size, std::vector<Int>* vec) const {
  if (!supernode_size) {
    return;
  }
  vec->reserve(vec->size() + supernode_size);
  Int index = principal_member;
  vec->push_back(index);
  while (index != tail_index_[principal_member]) {
    index = next_index_[index];
    vec->push_back(index);
  }
}

inline void QuotientGraph::ComputePostorder(std::vector<Int>* postorder) const {
  // Reconstruct the child links from the parent links in a contiguous array
  // (similar to the CSR format) by first counting the number of children of
  // each node.
  std::vector<Int> child_offsets(num_original_vertices_ + 1, 0);
  for (const Int& i : elimination_order_) {
    if (parents_[i] != -1) {
      ++child_offsets[parents_[i]];
    }
  }
  Int num_total_children = 0;
  for (Int i = 0; i <= num_original_vertices_; ++i) {
    const Int num_children = child_offsets[i];
    child_offsets[i] = num_total_children;
    num_total_children += num_children;
  }

  // Pack the children into a buffer.
  std::vector<Int> children(num_total_children);
  auto offsets_copy = child_offsets;
  for (const Int& i : elimination_order_) {
    const Int parent = parents_[i];
    if (parent != -1) {
      children[offsets_copy[parent]++] = i;
    }
  }

  // Scan for the roots and launch a pe-order traversal on each of them.
  postorder->resize(num_original_vertices_);
  std::vector<Int>::iterator iter = postorder->begin();
  for (const Int& index : elimination_order_) {
    if (!node_flags_[index]) {
      // This element was absorbed into another element.
      continue;
    }
    iter = PreorderTree(index, children, child_offsets, iter);
  }

  // Reverse the preordering (to form a postordering) in-place.
  std::reverse(postorder->begin(), iter);

  if (num_dense_ > 0) {
    // Pack the dense columns into the back.
    for (Int i = 0; i < num_original_vertices_; ++i) {
      if (!signed_supernode_sizes_[i] && parents_[i] == -1) {
        *iter = i;
        ++iter;
      }
    }
  }

  QUOTIENT_ASSERT(iter == postorder->end(),
      "Postorder had incorrect final offset.");
}

inline std::vector<Int>::iterator QuotientGraph::PreorderTree(
    Int root,
    const std::vector<Int>& children,
    const std::vector<Int>& child_offsets,
    std::vector<Int>::iterator iter) const {
  std::vector<Int> stack;
  stack.reserve(num_original_vertices_);

  stack.push_back(root);

  while (!stack.empty()) {
    // Pop a principal variable from the stack.
    const Int element = stack.back(); stack.pop_back();

    // Push the supernode into the preorder.
    {
      Int i = element;
      *(iter++) = element;
      while (i != tail_index_[element]) {
        i = next_index_[i];
        *(iter++) = i;
      }
    }

    // Push the children onto the stack. 
    for (Int index = child_offsets[element];
        index < child_offsets[element + 1]; ++index) {
      stack.push_back(children[index]);
    }
  }

  return iter;
}

inline const std::vector<Int>& QuotientGraph::Parents() const {
  return parents_;
}

inline void QuotientGraph::ExternalDegrees() {
  QUOTIENT_START_TIMER(timers_, kExternalDegrees);
  const Int shift = external_degree_shift_;
  const bool aggressive_absorption = control_.aggressive_absorption;

  for (const Int& i : elements_[pivot_]) {
    const Int supernode_size = -signed_supernode_sizes_[i];
    QUOTIENT_ASSERT(supernode_size > 0,
        "supernode " + std::to_string(i) + " had non-positive signed size "
        "when computing external element sizes");
    const Int shift_minus_supernode_size = shift - supernode_size;

    const Int offset = element_list_offsets_[i];
    const Int num_elements = element_list_sizes_[i];
    for (Int k = offset; k < offset + num_elements; ++k) {
      const Int element = element_and_adjacency_lists_[k];
      Int& shifted_external_degree = node_flags_[element];

      if (shifted_external_degree >= shift) {
        QUOTIENT_ASSERT(parents_[element] == -1,
          "Tried to subtract from an absorbed element.");
        QUOTIENT_ASSERT(shifted_external_degree != shift,
          "Shifted external size was equal to shift before subtracting");
        shifted_external_degree -= supernode_size;
      } else if (shifted_external_degree) {
        QUOTIENT_ASSERT(parents_[element] == -1,
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
        parents_[element] = pivot_;
        SwapClearVector(&elements_[element]);
        element_list_sizes_[element] = 0;
      }
    }
  }

  QUOTIENT_STOP_TIMER(timers_, kExternalDegrees);
}

inline void QuotientGraph::ResetExternalDegrees() {
  if (external_degree_shift_ + max_degree_ < max_shift_value_) {
    external_degree_shift_ += max_degree_ + 1;
    return; 
  }

  // Reset all non-absorbed shifted degrees to 1. Afterwards, the absorbed
  // elements will be marked as 0, the unabsorbed as 1, and the shift will be
  // 2.
#ifdef QUOTIENT_DEBUG
  std::cerr << "Resetting external degrees." << std::endl;
#endif
  external_degree_shift_ = 2;
  for (UInt i = 0; i < node_flags_.size(); ++i) {
    if (node_flags_[i]) {
      node_flags_[i] = 1;
    }
  }
}

inline std::vector<std::pair<std::string, double>>
QuotientGraph::ComponentSeconds() const {
  std::vector<std::pair<std::string, double>> times;
#ifdef QUOTIENT_ENABLE_TIMERS
  for (const std::pair<std::string, Timer>& pairing : timers_) {
    times.emplace_back(pairing.first, pairing.second.TotalSeconds());
  }
#endif
  return times;
}

} // namespace quotient

#endif // ifndef QUOTIENT_QUOTIENT_GRAPH_IMPL_H_
