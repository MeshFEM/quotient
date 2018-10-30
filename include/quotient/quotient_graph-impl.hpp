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

template<typename T>
void PrintVector(const std::vector<T>& vec, const std::string& msg) {
  std::cout << msg << ": ";
  for (std::size_t i = 0; i < vec.size(); ++i) {
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
  external_element_size_shift_(0),
  max_element_size_(0),
  max_shift_value_(std::numeric_limits<Int>::max() - graph.NumSources()),
  num_hash_bucket_collisions_(0),
  num_hash_collisions_(0),
  num_aggressive_absorptions_(0) {
  // Initialize the supernodes as simple.
  supernode_sizes_.resize(num_original_vertices_, 1);
  elimination_order_.reserve(num_original_vertices_);
  parents_.resize(num_original_vertices_, -1);
  next_index_.resize(num_original_vertices_, -1);
  tail_index_.resize(num_original_vertices_);
  std::iota(tail_index_.begin(), tail_index_.end(), 0);

  // Initialize the adjacency lists from the original graph.
  adjacency_lists_.resize(num_original_vertices_);
  const std::vector<GraphEdge>& edges = graph.Edges();
  for (Int source = 0; source < num_original_vertices_; ++source) {
    const Int source_edge_offset = graph.SourceEdgeOffset(source);
    const Int next_source_edge_offset = graph.SourceEdgeOffset(source + 1);
    const Int num_connections = next_source_edge_offset - source_edge_offset;
    std::vector<Int>& adjacency_list = adjacency_lists_[source];

    adjacency_list.reserve(num_connections);
    for (Int edge_index = source_edge_offset;
         edge_index < next_source_edge_offset; ++edge_index) {
      const GraphEdge& edge = edges[edge_index];
      if (edge.second != source) {
        adjacency_list.push_back(edge.second);
      }
    }
  }

  // Initialize the degree lists.
  degree_lists_.degrees.resize(num_original_vertices_);
  degree_lists_.heads.resize(num_original_vertices_ - 1, -1);
  degree_lists_.next_member.resize(num_original_vertices_, -1);
  degree_lists_.last_member.resize(num_original_vertices_, -1);
  for (Int source = 0; source < num_original_vertices_; ++source) {
    const Int degree = adjacency_lists_[source].size();
    degree_lists_.AddDegree(source, degree);
  }

  // Initialize the hash lists.
  hash_lists_.buckets.resize(num_original_vertices_);
  hash_lists_.hashes.resize(num_original_vertices_);
  hash_lists_.heads.resize(num_original_vertices_, -1);
  hash_lists_.next_member.resize(num_original_vertices_, -1);

  // Trivially initialize the element lists.
  element_lists_.resize(num_original_vertices_);

  // Trivially initialize the lower-triangular nonzero structures.
  if (control_.store_structures) {
    structures_.resize(num_original_vertices_);
  }
  elements_.resize(num_original_vertices_);
  element_sizes_.resize(num_original_vertices_, 0);

  shifted_external_element_sizes_.resize(num_original_vertices_, -1);

#ifdef QUOTIENT_DEBUG
  const bool using_exact_degree_mask = true;
#else
  const bool using_exact_degree_mask =
      control_.degree_type == kExactExternalDegree;
#endif
  if (using_exact_degree_mask) {
    exact_degree_mask_.resize(num_original_vertices_, 0);
  }
}

inline const std::vector<Int>& QuotientGraph::EliminationOrder() const {
  return elimination_order_;
}

inline Int QuotientGraph::GetNextPivot() {
  pivot_ = degree_lists_.FindMinimalIndex(control_.force_minimal_pivot_indices);
  return pivot_;
}

inline Int QuotientGraph::NumAggressiveAbsorptions() const {
  return num_aggressive_absorptions_;
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
    if (supernode_sizes_[i] == 0) {
      continue;
    }
    std::cout << "Supernode " << i << "\n";
    const std::vector<Int> supernode = FormSupernode(i);
    PrintVector(supernode, "  members");
    PrintVector(elements_[i], "  elements");
    PrintVector(adjacency_lists_[i], "  adjacency_list");
    PrintVector(element_lists_[i], "  element_list");
    std::cout << "\n";
  }
}

inline std::vector<Int> QuotientGraph::FormSupernode(
    Int principal_member) const {
  std::vector<Int> supernode;
  const Int supernode_size = std::abs(supernode_sizes_[principal_member]);
  AppendSupernode(principal_member, supernode_size, &supernode);
  return supernode;
}

inline const std::vector<Int>& QuotientGraph::Element(
    Int principal_member) const {
  return elements_[principal_member];
}

inline void QuotientGraph::FormEliminatedStructures(
    std::vector<std::vector<Int>>* eliminated_structures) const {
  eliminated_structures->resize(elimination_order_.size());
  for (std::size_t index = 0; index < elimination_order_.size(); ++index) {
    std::vector<Int>& eliminated_structure = (*eliminated_structures)[index];
    eliminated_structure = structures_[elimination_order_[index]];
  }
}

inline void QuotientGraph::ComputePivotStructure() {
  QUOTIENT_ASSERT(elements_[pivot_].empty(), "Chose a pivot more than once.");
  supernode_sizes_[pivot_] *= -1;

  // Push the supervariables in the pivot adjacency list into the structure.
  // TODO(Jack Poulson): Reserve space for elements_[pivot_].
  Int pivot_element_size = 0;
  for (const Int& index : adjacency_lists_[pivot_]) {
    const Int supernode_size = supernode_sizes_[index];
    QUOTIENT_ASSERT(supernode_size >= 0,
        "Encountered element in adjacency list.");
    if (!supernode_size) {
      continue;
    }
    QUOTIENT_ASSERT(parents_[index] == -1,
        "Absorbed element in adjacency list.");

    elements_[pivot_].push_back(index);
    pivot_element_size += supernode_size;
    supernode_sizes_[index] = -supernode_size;

    if (control_.store_structures) {
      AppendSupernode(index, supernode_size, &structures_[pivot_]);
    }
  }

  // Push the unique supervariables in the patterns of the pivot element list
  // into the current structure.
  Int num_kept_elements = 0;
  std::vector<Int>& pivot_elements = element_lists_[pivot_];
  for (std::size_t elem_index = 0; elem_index < pivot_elements.size();
      ++elem_index) {
    const Int element = pivot_elements[elem_index];
    QUOTIENT_ASSERT(parents_[element] == -1,
        "Used an absorbed element in pivot structure.");
    parents_[element] = pivot_;
    pivot_elements[num_kept_elements++] = element;

    for (const Int& index : elements_[element]) {
      const Int supernode_size = supernode_sizes_[index];
      // While no eliminated supernodes should appear in unabsorbed elements,
      // we have (temporarily) flipped the signs of the supernode sizes of
      // the members we have already added to this pivot's structure.
      if (supernode_size <= 0) {
        continue;
      }

      elements_[pivot_].push_back(index);
      pivot_element_size += supernode_size;
      supernode_sizes_[index] = -supernode_size;

      if (control_.store_structures) {
        AppendSupernode(index, supernode_size, &structures_[pivot_]);
      }
    }
  }
  pivot_elements.resize(num_kept_elements);

  element_sizes_[pivot_] = pivot_element_size;
  if (pivot_element_size > max_element_size_) {
    max_element_size_ = pivot_element_size;
  }

#ifdef QUOTIENT_DEBUG
  for (const Int& element : element_lists_[pivot_]) {
    QUOTIENT_ASSERT(parents_[element] == pivot_,
        "Pivot element's parent was not properly marked.");

    Int element_size = 0;
    for (const Int& index : elements_[element]) {
      const Int supernode_size = -supernode_sizes_[index];
      QUOTIENT_ASSERT(supernode_size >= 0,
          "Flipped supernode size was expected to be positive.");
      element_size += supernode_size;
    }
    QUOTIENT_ASSERT(element_size == element_sizes_[element],
      "Element size did not match its cached value.");
  }
#endif
}

inline Int QuotientGraph::NumPivotElements() const {
  return element_lists_[pivot_].size();
}

inline Int QuotientGraph::NumPivotDegreeUpdates() const {
  return elements_[pivot_].size();
}

inline Int QuotientGraph::NumPivotDegreeUpdatesWithMultipleElements() const {
  Int num_multi_updates = 0;
  for (const Int& i : elements_[pivot_]) {
    if (element_lists_[i].size() > 2) {
      ++num_multi_updates;
    }
  }
  return num_multi_updates;
}

inline Int QuotientGraph::NumPivotCholeskyNonzeros() const {
  const Int pivot_size = std::abs(supernode_sizes_[pivot_]);
  const Int structure_size = element_sizes_[pivot_];
  const Int diag_block_nonzeros = (pivot_size * (pivot_size + 1)) / 2;
  const Int subdiagonal_nonzeros = structure_size * pivot_size;
  return diag_block_nonzeros + subdiagonal_nonzeros;
}

inline double QuotientGraph::NumPivotCholeskyFlops() const {
  const Int pivot_size = std::abs(supernode_sizes_[pivot_]);
  const Int structure_size = element_sizes_[pivot_];
  const double diag_block_flops =  std::pow(1. * pivot_size, 3.) / 3.;
  const double schur_complement_flops =
      std::pow(1. * structure_size, 2.) * pivot_size;
  return diag_block_flops + schur_complement_flops;
}

inline std::pair<Int, std::size_t> QuotientGraph::PackCountAndHashAdjacencies(
    Int principal_variable) {
  Int degree = 0;
  std::size_t hash = 0;

  Int num_packed = 0;
  std::vector<Int>& adjacency_list = adjacency_lists_[principal_variable];
  const Int num_adjacencies = adjacency_list.size();
  for (Int index = 0; index < num_adjacencies; ++index) {
    const Int j = adjacency_list[index];
    // We must filter out: non-principal variables, eliminated supernodes,
    // the pivot supernode, and all members of the pivot structure.
    //
    // Due to our (temporary) negation of the supernode sizes of the pivot and
    // members of the pivot structure, all four cases can be handled by
    // demanding a positive (signed) supernode size.
    const Int supernode_size = supernode_sizes_[j];
    if (supernode_size > 0) {
      degree += supernode_size;
      QUOTIENT_HASH_COMBINE(hash, j);
      adjacency_list[num_packed++] = j;
    }
  }
  adjacency_list.resize(num_packed);

  return std::make_pair(degree, hash);
}

inline std::pair<Int, std::size_t>
QuotientGraph::ExactEmptyExternalDegreeAndHash(Int principal_variable) {
  QUOTIENT_ASSERT(element_lists_[principal_variable].empty(),
      "The element list was falsely assumed empty.");

  // Add the cardinality of A_i \ supernode(i), where 'i' is the principal
  // variable.
  std::pair<Int, std::size_t> degree_and_hash =
      PackCountAndHashAdjacencies(principal_variable);

  return degree_and_hash;
}

inline std::pair<Int, std::size_t>
QuotientGraph::ExactSingleExternalDegreeAndHash(Int principal_variable) {
  // We should only have one member of the element list, 'pivot'.
  QUOTIENT_ASSERT(element_lists_[principal_variable].size() == 1,
      "The element list was falsely assumed to have a single entry.");
  QUOTIENT_ASSERT(element_lists_[principal_variable][0] == pivot_,
      "The element list should have only contained the pivot.");

  // Add the cardinality of A_i \ supernode(i), where 'i' is the principal
  // variable.
  std::pair<Int, std::size_t> degree_and_hash =
      PackCountAndHashAdjacencies(principal_variable);

  // Add |L_p \ supernode(i)|.
  const Int supernode_size = -supernode_sizes_[principal_variable];
  QUOTIENT_ASSERT(supernode_size > 0,
      "The flipped supernode size should have been positive.");
  degree_and_hash.first += element_sizes_[pivot_] - supernode_size;
  QUOTIENT_HASH_COMBINE(degree_and_hash.second, pivot_);

  return degree_and_hash;
}

inline std::pair<Int, std::size_t>
QuotientGraph::ExactDoubleExternalDegreeAndHash(Int principal_variable) {
  const bool aggressive_absorption = control_.aggressive_absorption;
  std::vector<Int>& element_list = element_lists_[principal_variable];

  // There should be exactly two members in the element list, and the second
  // should be the pivot.
  QUOTIENT_ASSERT(element_list.size() == 2,
      "The element list should have had exactly two members.");
  QUOTIENT_ASSERT(element_list[1] == pivot_,
      "The element list should have contained the pivot as its 2nd member.");

  // Add the cardinality of A_i \ supernode(i), where 'i' is the principal
  // variable.
  std::pair<Int, std::size_t> degree_and_hash =
      PackCountAndHashAdjacencies(principal_variable);

  // Add |L_p \ supernode(i)|.
  const Int supernode_size = -supernode_sizes_[principal_variable];
  QUOTIENT_ASSERT(supernode_size > 0,
      "The flipped supernode size should have been positive.");
  degree_and_hash.first += element_sizes_[pivot_] - supernode_size;
  QUOTIENT_HASH_COMBINE(degree_and_hash.second, pivot_);

  const Int element = element_list[0];
  const Int external_element_size =
      shifted_external_element_sizes_[element] - external_element_size_shift_;
  if (!aggressive_absorption || external_element_size != 0) {
    // Add |L_e \ L_p|.
    degree_and_hash.first += external_element_size;
    QUOTIENT_HASH_COMBINE(degree_and_hash.second, element);
  } else {
    element_list[0] = pivot_;
    element_list.pop_back();
  }

  return degree_and_hash;
}

inline std::pair<Int, std::size_t>
QuotientGraph::ExactGenericExternalDegreeAndHash(Int principal_variable) {
  const bool aggressive_absorption = control_.aggressive_absorption;
  std::vector<Int>& element_list = element_lists_[principal_variable];

  // Add the cardinality of A_i \ supernode(i), where 'i' is the principal
  // variable.
  std::pair<Int, std::size_t> degree_and_hash =
      PackCountAndHashAdjacencies(principal_variable);

  // Add on the number of unique entries in the structures of the element lists
  // that are outside supernode(i).
  Int num_packed = 0;
  for (const Int& element : element_list) {
    if (aggressive_absorption) {
      if (parents_[element] != -1) {
        // Skip absorbed elements. This should only be possible with
        // aggressive absorption (so the external element size should be 0).
        continue;
      }
      element_list[num_packed++] = element;
    }
    QUOTIENT_HASH_COMBINE(degree_and_hash.second, element);

    for (const Int& j : elements_[element]) {
      // Unabsorbed elements should not have any eliminated members of
      // in their element. Thus, we can take the absolute value of the signed
      // supernode size (since some members might be in the pivot structure).
      if (exact_degree_mask_[j] || principal_variable == j) {
        continue;
      }
      degree_and_hash.first += std::abs(supernode_sizes_[j]);
      exact_degree_mask_[j] = 1;
    }
  }
  if (aggressive_absorption) {
    element_list.resize(num_packed);
  }

  // Clear the mask.
  // NOTE: We allow (*mask)[i] to be set to zero to avoid branching.
  for (const Int& element : element_list) {
    for (const Int& j : elements_[element]) {
      exact_degree_mask_[j] = 0;
    }
  }

  return degree_and_hash;
}

inline std::pair<Int, std::size_t>
QuotientGraph::ExactExternalDegreeAndHash(Int principal_variable) {
  const Int num_elements = element_lists_[principal_variable].size();
  if (num_elements == 0) {
    return ExactEmptyExternalDegreeAndHash(principal_variable);
  } else if (num_elements == 1) {
    return ExactSingleExternalDegreeAndHash(principal_variable);
  } else if (num_elements == 2) {
    return ExactDoubleExternalDegreeAndHash(principal_variable);
  } else {
    return ExactGenericExternalDegreeAndHash(principal_variable);
  }
}

inline std::pair<Int, std::size_t>
QuotientGraph::AmestoyExternalDegreeAndHash(Int principal_variable) {
  const bool aggressive_absorption = control_.aggressive_absorption;
  const Int shift = external_element_size_shift_;
  const Int supernode_size = -supernode_sizes_[principal_variable];
  QUOTIENT_ASSERT(supernode_size > 0,
      "The negated supernode size was expected to be positive.");
  std::vector<Int>& element_list = element_lists_[principal_variable];

  // Note that this usage of 'external' refers to |L_p \ supernode(i)| and not
  // |L_e \ L_p|, as is the case for 'external_element_sizes', where 'i'
  // is 'principal_variable'.
  const Int external_pivot_element_size =
      element_sizes_[pivot_] - supernode_size;
  QUOTIENT_ASSERT(external_pivot_element_size >= 0,
      "Encountered a negative external_pivot_element_size");

  const Int num_vertices_left =
      num_original_vertices_ - num_eliminated_vertices_;
  const Int bound0 = num_vertices_left - supernode_size;

  const Int old_degree = degree_lists_.degrees[principal_variable];
  const Int bound1 = old_degree + external_pivot_element_size;

  // bound_2 = |A_i \ supernode(i)| + |L_p \ supernode(i)| +
  //           \sum_{e in E_i \ {p}} |L_e \ L_p|.
  std::pair<Int, std::size_t> degree_and_hash =
      PackCountAndHashAdjacencies(principal_variable);
  degree_and_hash.first += external_pivot_element_size;

  Int num_packed = 0;
  const Int num_elements = element_list.size();
  QUOTIENT_ASSERT(element_list.back() == pivot_,
      "Pivot was not at the back of the element list.");
  for (Int elem_index = 0; elem_index < num_elements - 1; ++elem_index) {
    const Int element = element_list[elem_index];
    QUOTIENT_ASSERT(element != pivot_, "Iterated over pivot element.");
    if (aggressive_absorption) {
      if (parents_[element] != -1) {
        // Skip absorbed elements. This should only be possible with
        // aggressive absorption.
        continue;
      }
      element_list[num_packed++] = element;
    }
    const Int external_element_size =
        shifted_external_element_sizes_[element] - shift;
    QUOTIENT_ASSERT(external_element_size >= 0,
        "Ran into a missing entry in external_element_sizes_");
    degree_and_hash.first += external_element_size;
    QUOTIENT_HASH_COMBINE(degree_and_hash.second, element);
  }
  if (aggressive_absorption) {
    element_list.resize(num_packed);
    element_list.push_back(pivot_);
  }
  QUOTIENT_HASH_COMBINE(degree_and_hash.second, pivot_);

  if (bound0 < degree_and_hash.first) {
    degree_and_hash.first = bound0;
  }
  if (bound1 < degree_and_hash.first) {
    degree_and_hash.first = bound1;
  }
 
  return degree_and_hash;
}

inline std::pair<Int, std::size_t>
QuotientGraph::AshcraftExternalDegreeAndHash(Int principal_variable) {
  if (element_lists_[principal_variable].size() == 2) {
    return ExactDoubleExternalDegreeAndHash(principal_variable);
  }
  return GilbertExternalDegreeAndHash(principal_variable);
}

inline std::pair<Int, std::size_t>
QuotientGraph::GilbertExternalDegreeAndHash(Int principal_variable) {
  const bool aggressive_absorption = control_.aggressive_absorption;
  std::vector<Int>& element_list = element_lists_[principal_variable];

  std::pair<Int, std::size_t> degree_and_hash =
      PackCountAndHashAdjacencies(principal_variable);

  const Int supernode_size = -supernode_sizes_[principal_variable];
  QUOTIENT_ASSERT(supernode_size > 0,
      "The negated supernode size was expected to be positive.");

  Int num_packed = 0;
  for (const Int& element : element_list) {
    if (aggressive_absorption) {
      if (parents_[element] != -1) {
        // Skip absorbed elements. This should only be possible with
        // aggressive absorption.
        continue;
      }
      element_list[num_packed++] = element;
    }
    QUOTIENT_ASSERT(parents_[element] == -1,
        "Used absorbed element in Gilbert degree update.");
    QUOTIENT_ASSERT(element_sizes_[element] - supernode_size >= 0,
        "Negative Gilbert degree update.");
    degree_and_hash.first += element_sizes_[element] - supernode_size;
    QUOTIENT_HASH_COMBINE(degree_and_hash.second, element);
  }
  if (aggressive_absorption) {
    element_list.resize(num_packed);
  }

  const Int num_vertices_left =
      num_original_vertices_ - num_eliminated_vertices_;
  const Int bound = num_vertices_left - supernode_size;
  if (bound < degree_and_hash.first) {
    degree_and_hash.first = bound;
  }

  return degree_and_hash;
}

inline std::pair<Int, std::size_t>
QuotientGraph::ExternalDegreeAndHash(Int principal_variable) {
  std::pair<Int, std::size_t> degree_and_hash;
  switch(control_.degree_type) {
    case kExactExternalDegree: {
      degree_and_hash = ExactExternalDegreeAndHash(principal_variable);
      break;
    }
    case kAmestoyExternalDegree: {
      degree_and_hash = AmestoyExternalDegreeAndHash(principal_variable);
      break;
    }
    case kAshcraftExternalDegree: {
      degree_and_hash = AshcraftExternalDegreeAndHash(principal_variable);
      break;
    }
    case kGilbertExternalDegree: {
      degree_and_hash = GilbertExternalDegreeAndHash(principal_variable);
      break;
    }
  }
  return degree_and_hash;
}

inline void QuotientGraph::ComputeExternalDegreesAndHashes(
    std::vector<Int>* external_degrees, std::vector<std::size_t>* bucket_keys) {
  const std::vector<Int>& pivot_element = elements_[pivot_];
  const std::size_t supernodal_struct_size = pivot_element.size();
  external_degrees->resize(supernodal_struct_size);
  bucket_keys->resize(supernodal_struct_size);
  // TODO(Jack Poulson): Consider how to parallelize using different masks
  // for each thread.
  for (std::size_t index = 0; index < supernodal_struct_size; ++index) {
    const Int i = pivot_element[index];
#ifdef QUOTIENT_DEBUG
    std::vector<std::pair<Int, std::size_t>> degrees_and_hashes(4);
    ExternalDegreeType degree_type = control_.degree_type;
    for (Int type = 0; type < 4; ++type) {
      control_.degree_type = static_cast<ExternalDegreeType>(type);
      degrees_and_hashes[type] = ExternalDegreeAndHash(i);
    }
    control_.degree_type = degree_type;
    QUOTIENT_ASSERT(
        degrees_and_hashes[1].first >= degrees_and_hashes[0].first &&
        degrees_and_hashes[2].first >= degrees_and_hashes[1].first &&
        degrees_and_hashes[3].first >= degrees_and_hashes[2].first,
        "Degrees (" + std::to_string(element_lists_[i].size()) + "): " +
        std::to_string(degrees_and_hashes[0].first) + ", " +
        std::to_string(degrees_and_hashes[1].first) + ", " +
        std::to_string(degrees_and_hashes[2].first) + ", " +
        std::to_string(degrees_and_hashes[3].first));
    (*external_degrees)[index] = degrees_and_hashes[degree_type].first;
    (*bucket_keys)[index] = degrees_and_hashes[degree_type].second;
#else
    const std::pair<Int, std::size_t> degree_and_hash =
        ExternalDegreeAndHash(i);
    (*external_degrees)[index] = degree_and_hash.first;
    (*bucket_keys)[index] = degree_and_hash.second;
#endif
  }
}

inline void QuotientGraph::UpdateExternalDegrees(
    const std::vector<Int>& external_degrees) {
  const std::vector<Int>& pivot_element = elements_[pivot_];
  const Int supernodal_struct_size = pivot_element.size();
  for (Int index = 0; index < supernodal_struct_size; ++index) {
    const Int i = pivot_element[index];
    const Int degree = external_degrees[index];
    degree_lists_.UpdateDegree(i, degree);
    QUOTIENT_ASSERT(degree_lists_.degrees[i] == degree,
        "Degree was not updated.");
  }
}

inline bool QuotientGraph::StructuralSupervariablesAreQuotientIndistinguishable(
    Int i, Int j) const {
  QUOTIENT_ASSERT(i != j,
      "Explicitly tested for equivalence of a supernode with itself...");

#ifdef QUOTIENT_DEBUG
  for (const Int& element : element_lists_[i]) {
    QUOTIENT_ASSERT(parents_[element] == -1,
        "Absorbed element was in element list during absorption test.");
  }
  for (const Int& element : element_lists_[j]) {
    QUOTIENT_ASSERT(parents_[element] == -1,
        "Absorbed element was in element list during absorption test.");
  }
  for (const Int& index : adjacency_lists_[i]) {
    QUOTIENT_ASSERT(supernode_sizes_[index] > 0,
        "Non-positive supernode size in adj list during absorption test.");
  }
  for (const Int& index : adjacency_lists_[j]) {
    QUOTIENT_ASSERT(supernode_sizes_[index] > 0,
        "Non-positive supernode size in adj list during absorption test.");
  }
#endif

  // Check if A_i = A_j.
  if (adjacency_lists_[i].size() != adjacency_lists_[j].size()) {
    return false;
  }
  for (std::size_t index = 0; index < adjacency_lists_[i].size(); ++index) {
    if (adjacency_lists_[i][index] != adjacency_lists_[j][index]) {
      return false;
    }
  }

  // Check if E_i = E_j.
  if (element_lists_[i].size() != element_lists_[j].size()) {
    return false;
  }
  // NOTE: The element lists are *not* sorted, but, if they are the same, they
  // should appear in the same order.
  for (std::size_t index = 0; index < element_lists_[i].size(); ++index) {
    if (element_lists_[i][index] != element_lists_[j][index]) {
#ifdef QUOTIENT_DEBUG
      std::vector<Int> elem_list_i = element_lists_[i];
      std::vector<Int> elem_list_j = element_lists_[j];
      std::sort(elem_list_i.begin(), elem_list_i.end());
      std::sort(elem_list_j.begin(), elem_list_j.end());
      for (std::size_t k = 0; k < elem_list_i.size(); ++k) {
        if (elem_list_i[k] != elem_list_j[k]) {
          return false;
        }
      }
      std::cerr << "Skipped an unsorted exact match." << std::endl;
#endif
      return false;
    }
  }

  return true;
}

inline void QuotientGraph::MergeVariables(
    const std::vector<std::size_t>& bucket_keys) {
  const std::vector<Int>& pivot_element = elements_[pivot_];
  const Int supernodal_struct_size = pivot_element.size();

  // Add the hashes into the hash lists.
  for (Int i_index = 0; i_index < supernodal_struct_size; ++i_index) {
    const std::size_t hash = bucket_keys[i_index];
    const Int bucket = hash % num_original_vertices_;
    const Int i = pivot_element[i_index];
    hash_lists_.AddHash(i, hash, bucket);
  }

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
      if (!supernode_sizes_[i]) {
        continue;
      }
      QUOTIENT_ASSERT(supernode_sizes_[i] < 0,
          "Supernode size should have been temporarily negative.");
      const std::size_t i_hash = hash_lists_.hashes[i];
      for (Int j = next_member[i]; j != -1; j = next_member[j]) {
        if (!supernode_sizes_[j]) { 
          continue;
        }
        QUOTIENT_ASSERT(supernode_sizes_[j] < 0,
            "Supernode size should have been temporarily negative.");
        const std::size_t j_hash = hash_lists_.hashes[j];
        if (i_hash != j_hash) {
          ++num_hash_bucket_collisions_;
          continue;
        }

        if (StructuralSupervariablesAreQuotientIndistinguishable(i, j)) {
          const Int absorbed_size = -supernode_sizes_[j];
          QUOTIENT_ASSERT(absorbed_size > 0,
              "Absorbed size should have been positive.");
 
          //#pragma omp critical
          const Int old_degree = degree_lists_.degrees[i];
          const Int degree = old_degree - absorbed_size;
          degree_lists_.UpdateDegree(i, degree);
          degree_lists_.RemoveDegree(j);

          // Merge [i] -> [j] (and i becomes the principal member).
          next_index_[tail_index_[i]] = j;
          tail_index_[i] = tail_index_[j];
          supernode_sizes_[i] -= absorbed_size;  // Recall it is negated.
          supernode_sizes_[j] = 0;

          SwapClearVector(&adjacency_lists_[j]);
          SwapClearVector(&element_lists_[j]);
        } else {
          ++num_hash_bucket_collisions_;
          ++num_hash_collisions_;
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
}

inline void QuotientGraph::ConvertPivotIntoElement(
    const std::vector<Int>& aggressive_absorption_elements) {
  const Int supernode_size = -supernode_sizes_[pivot_];
  QUOTIENT_ASSERT(supernode_size > 0,
      "The supernode size was assumed positive.");
  const std::vector<Int>& pivot_element_list = element_lists_[pivot_];

  // Since this supervariable is being eliminated, it needs to be implicitly
  // removed from all elements containing it. This is accomplished through
  // decrementing their element sizes and marking the supernode as eliminated
  // by making its sign negative.
  for (const Int& element : pivot_element_list) {
    element_sizes_[element] -= supernode_size;
  }

  // Unflip the signs of supernode_sizes for the members of the pivot structure
  // and simultaneously remove the non-principal members.
  Int num_packed = 0;
  for (const Int& i : elements_[pivot_]) {
    QUOTIENT_ASSERT(supernode_sizes_[i] <= 0,
        "A member of pivot element had a positive signed supernode size.");
    if (supernode_sizes_[i] < 0) {
      supernode_sizes_[i] = -supernode_sizes_[i];
      elements_[pivot_][num_packed++] = i;
    }
  }
  elements_[pivot_].resize(num_packed);

  degree_lists_.RemoveDegree(pivot_);
  SwapClearVector(&adjacency_lists_[pivot_]);

  // Store the bidirectional links for element absorptions so that we can later
  // quickly compute the post-ordering.
  if (!aggressive_absorption_elements.empty()) {
    // Append the aggressive absorption elements onto the element list so that
    // the element list becomes the list of *absorbed* elements.
    std::vector<Int>& element_list = element_lists_[pivot_];
    element_list.reserve(
        element_list.size() + aggressive_absorption_elements.size());
    for (const Int& element : aggressive_absorption_elements) {
#ifdef QUOTIENT_DEBUG
      for (const Int& e : element_list) {
        QUOTIENT_ASSERT(e != element,
            "Aggressive absorption element was redundant.");
      }
#endif
      element_list.push_back(element); 
    }
  }

  elimination_order_.push_back(pivot_);
  num_eliminated_vertices_ += supernode_size;
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

inline void QuotientGraph::ComputePreorder(std::vector<Int>* preorder) const {
  // Scan for the roots and launch a pre-order traversal on each of them.
  preorder->resize(num_original_vertices_);
  std::vector<Int>::iterator iter = preorder->begin();
  for (const Int& index : elimination_order_) {
    if (parents_[index] != -1) {
      continue;
    }
    iter = PreorderTree(index, iter);
  }
  QUOTIENT_ASSERT(iter == preorder->end(),
      "Preorder had incorrect final offset.");
}

inline std::vector<Int>::iterator QuotientGraph::PreorderTree(
    Int root, std::vector<Int>::iterator iter) const {
  std::vector<Int> stack;
  stack.reserve(num_original_vertices_);

  stack.push_back(root);

  while (!stack.empty()) {
    // Pop a principal variable from the stack.
    const Int element = stack.back(); stack.pop_back();

    // Push the supernode into the preorder.
    {
      Int index = element;
      *(iter++) = element;
      while (index != tail_index_[element]) {
        index = next_index_[index];
        *(iter++) = index;
      }
    }

    // Push the children onto the stack. 
    for (const Int& e : element_lists_[element]) {
      stack.push_back(e); 
    }
  }

  return iter;
}

inline void QuotientGraph::AbsorptionAndExternalElementSizes(
    std::vector<Int>* aggressive_absorption_elements) {
  const Int shift = external_element_size_shift_;
  const bool aggressive_absorption = control_.aggressive_absorption;
  if (aggressive_absorption) {
    // Follow the advice at the beginning of Section 5 of [AMD-96] and absorb
    // any element e that satisfies |L_e \ L_p| = 0.
    aggressive_absorption_elements->clear();
  }

  for (const Int& i : elements_[pivot_]) {
    std::vector<Int>& element_list = element_lists_[i];
    const Int supernode_i_size = -supernode_sizes_[i];
    QUOTIENT_ASSERT(supernode_i_size > 0,
        "supernode " + std::to_string(i) +
        " had non-positive signed size when computing external element sizes");

    Int num_packed = 0;
    const Int num_elements = element_list.size();
    for (Int elem_index = 0; elem_index < num_elements; ++elem_index) {
      const Int element = element_list[elem_index];
      if (parents_[element] != -1) {
        continue;
      }

      Int& shifted_external_size = shifted_external_element_sizes_[element];
      if (shifted_external_size < shift) {
        shifted_external_size =
            (element_sizes_[element] - supernode_i_size) + shift;
      } else {
        shifted_external_size -= supernode_i_size;
      }
      QUOTIENT_ASSERT(
          shifted_external_size >= shift,
          "Computed negative external element size.");

      // Mark any element with no exterior element size that is not a member
      // of the pivot element list for absorption if aggressive absorption is
      // enabled.
      if (aggressive_absorption && shifted_external_size == shift) {
        ++num_aggressive_absorptions_;
        parents_[element] = pivot_;
        aggressive_absorption_elements->push_back(element);
      } else {
        element_list[num_packed++] = element;
      }
    }
    element_list.resize(num_packed);
    element_list.push_back(pivot_);
  }
}

inline void QuotientGraph::ResetExternalElementSizes() {
  if (external_element_size_shift_ + max_element_size_ < max_shift_value_) {
    external_element_size_shift_ += max_element_size_ + 1;
    return; 
  }

  for (std::size_t i = 0; i < shifted_external_element_sizes_.size(); ++i) {
    shifted_external_element_sizes_[i] = -1;
  }
  external_element_size_shift_ = 0;
}

} // namespace quotient

#endif // ifndef QUOTIENT_QUOTIENT_GRAPH_IMPL_H_
