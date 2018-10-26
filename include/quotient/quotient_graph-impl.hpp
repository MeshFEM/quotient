/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_QUOTIENT_GRAPH_IMPL_H_
#define QUOTIENT_QUOTIENT_GRAPH_IMPL_H_

#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>

#include "quotient/config.hpp"
#include "quotient/coordinate_graph.hpp"
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
  num_hash_collisions_(0) {
  // Initialize the supernodes as simple.
  supernode_sizes_.resize(num_original_vertices_, 1);
  elimination_order_.reserve(num_original_vertices_);
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
  degree_lists_.degree_heads.resize(num_original_vertices_ - 1, -1);
  degree_lists_.next_degree_member.resize(num_original_vertices_, -1);
  degree_lists_.last_degree_member.resize(num_original_vertices_, -1);
  for (Int source = 0; source < num_original_vertices_; ++source) {
    const Int degree = adjacency_lists_[source].size();
    degree_lists_.AddDegree(source, degree);
  }

  // Trivially initialize the element lists.
  element_lists_.resize(num_original_vertices_);

  // Trivially initialize the lower-triangular nonzero structures.
  if (control_.store_structures) {
    structures_.resize(num_original_vertices_);
  }
  elements_.resize(num_original_vertices_);
  element_sizes_.resize(num_original_vertices_, 0);

  buckets_.resize(num_original_vertices_);

#ifdef QUOTIENT_DEBUG
  // We will compute all degree approximations to ensure that the guaranteed
  // inequalities hold.
  using_external_element_sizes_ = true;
#else
  using_external_element_sizes_ =
      control_.aggressive_absorption ||
      control_.degree_type == kExactExternalDegree ||
      control_.degree_type == kAmestoyExternalDegree ||
      control_.degree_type == kAshcraftExternalDegree;
#endif
  if (using_external_element_sizes_) {
    external_element_sizes_.resize(num_original_vertices_, -1);
  }

  pivot_mask_.resize(num_original_vertices_, 0);

#ifdef QUOTIENT_DEBUG
  exact_degree_mask_.resize(num_original_vertices_, 0);
#else
  if (control_.degree_type == kExactExternalDegree) {
    exact_degree_mask_.resize(num_original_vertices_, 0);
  }
#endif
}

inline const std::vector<Int>& QuotientGraph::EliminationOrder() const {
  return elimination_order_;
}

inline Int QuotientGraph::GetNextPivot() {
  pivot_ = degree_lists_.FindMinimalIndex(control_.force_minimal_pivot_indices);
  return pivot_;
}

inline bool QuotientGraph::UsingExternalElementSizes() const {
  return using_external_element_sizes_;
}

inline const std::vector<std::pair<Int, Int>>&
QuotientGraph::VariableMerges() const {
  return variable_merges_;
}

inline const std::vector<std::pair<Int, Int>>&
QuotientGraph::AggressiveAbsorptions() const {
  return aggressive_absorptions_;
}

inline Int QuotientGraph::NumOriginalVertices() const {
  return num_original_vertices_;
}

inline Int QuotientGraph::NumEliminatedVertices() const {
  return num_eliminated_vertices_;
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
    if (supernode_sizes_[i] < 0) {
      PrintVector(elements_[i], "  elements");
    } else {
      PrintVector(adjacency_lists_[i], "  adjacency_list");
      PrintVector(element_lists_[i], "  element_list");
    }
    std::cout << "\n";
  }
}

inline std::vector<Int> QuotientGraph::FormSupernode(Int i) const {
  std::vector<Int> supernode;
  if (supernode_sizes_[i] == 0) {
    return supernode;
  }

  supernode.reserve(std::abs(supernode_sizes_[i]));
  Int index = i;
  supernode.push_back(index);
  while (index != tail_index_[i]) {
    index = next_index_[index];
    supernode.push_back(index);
  }

  return supernode;
}

inline void QuotientGraph::FormEliminatedStructures(
    std::vector<std::vector<Int>>* eliminated_structures) const {
  eliminated_structures->resize(elimination_order_.size());
  for (std::size_t index = 0; index < elimination_order_.size(); ++index) {
    std::vector<Int>& eliminated_structure = (*eliminated_structures)[index];
    eliminated_structure = structures_[elimination_order_[index]];
    std::sort(eliminated_structure.begin(), eliminated_structure.end());
  }
}

inline Int QuotientGraph::ComputePivotStructure() {
#ifdef QUOTIENT_DEBUG
  for (Int index = 0; index < num_original_vertices_; ++index) {
    if (pivot_mask_[index]) {
      std::cerr << "Pivot mask had a nonzero on entry to ComputePivotStructure."
                << std::endl;
    }
  }
  if (!elements_[pivot_].empty()) {
    std::cerr << "Chose a pivot more than once." << std::endl;
  }
#endif
  // TODO(Jack Poulson): Reserve space for elements_[pivot_].
  for (const Int& index : adjacency_lists_[pivot_]) {
    const Int supernode_size = supernode_sizes_[index];
#ifdef QUOTIENT_DEBUG
    if (supernode_size < 0) {
      std::cerr << "Encountered an element in the adjacency list." << std::endl;
    }
#endif
    if (!supernode_size) {
      continue;
    }
    pivot_mask_[index] = 1;
    elements_[pivot_].push_back(index);
    element_sizes_[pivot_] += supernode_size;

    if (control_.store_structures) {
      // Fill in the supernode.
      Int k = index;
      structures_[pivot_].push_back(k);
      while (k != tail_index_[index]) {
        k = next_index_[k];
        structures_[pivot_].push_back(k);
      }
    }
  }

  Int num_stale_element_members = 0;
  for (const Int& element : element_lists_[pivot_]) {
    for (const Int& index : elements_[element]) {
      const Int supernode_size = supernode_sizes_[index];
      if (supernode_size < 0) {
        ++num_stale_element_members;
        // TODO(Jack Poulson): Quickly remove elements[element][index]. Though
        // this branch appears to be extremely rare.
      }
      if (index == pivot_ || supernode_size <= 0 || pivot_mask_[index]) {
        continue;
      }
      pivot_mask_[index] = 1;
      elements_[pivot_].push_back(index);
      element_sizes_[pivot_] += supernode_size;

      if (control_.store_structures) {
        // Fill in the supernode.
        Int k = index;
        structures_[pivot_].push_back(k);
        while (k != tail_index_[index]) {
          k = next_index_[k];
          structures_[pivot_].push_back(k);
        }
      }
    }
  }

  // If aggressive absorption is to be activated, we need to save a copy of
  // the pivot element list so that element_sizes can be later decremented after
  // converting the pivot into an element.
  if (control_.aggressive_absorption) {
    original_pivot_element_list_ = element_lists_[pivot_];
  }

#ifdef QUOTIENT_DEBUG
  for (const Int& element : element_lists_[pivot_]) {
    Int element_size = 0;
    for (const Int& index : elements_[element]) {
      const Int supernode_size = supernode_sizes_[index];
      if (supernode_size <= 0) {
        continue;
      }
      element_size += supernode_size;
    }
    const Int cached_element_size = element_sizes_[element];
    if (element_size != cached_element_size) {
      std::cerr << "Element " << element << " had a size of "
                << element_size << " but the cached size was "
                << cached_element_size << std::endl;
    }
  }
#endif

  return num_stale_element_members;
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
  const Int pivot_size = supernode_sizes_[pivot_];
  const Int structure_size = element_sizes_[pivot_];
  const Int diag_block_nonzeros = (pivot_size * (pivot_size + 1)) / 2;
  const Int subdiagonal_nonzeros = structure_size * pivot_size;
  return diag_block_nonzeros + subdiagonal_nonzeros;
}

inline double QuotientGraph::NumPivotCholeskyFlops() const {
  const Int pivot_size = supernode_sizes_[pivot_];
  const Int structure_size = element_sizes_[pivot_];
  const double diag_block_flops =  std::pow(1. * pivot_size, 3.) / 3.;
  const double schur_complement_flops =
      std::pow(1. * structure_size, 2.) * pivot_size;
  return diag_block_flops + schur_complement_flops;
}

inline void QuotientGraph::UpdateAdjacencyListsAfterSelectingPivot() {
  for (const Int& i : elements_[pivot_]) {
    // Remove redundant adjacency entries:
    //   A_i := (A_i \ L_p) \ supernode(p).
    Int packed_size = 0;
    for (Int index : adjacency_lists_[i]) {
      if (index == pivot_ || pivot_mask_[index] || !supernode_sizes_[index]) {
        continue;
      }
      adjacency_lists_[i][packed_size++] = index;
    }
    adjacency_lists_[i].resize(packed_size);
  }
}

inline void QuotientGraph::UpdateElementListsAfterSelectingPivot() {
  for (const Int& i : elements_[pivot_]) {
    // Element absorption:
    //   E_i := (E_i \ E_p) \cup {p}.
    std::vector<Int>& element_list = element_lists_[i];
    Int num_packed = 0;
    for (const Int& element : element_list) {
      if (pivot_mask_[element]) {
        continue;
      }
      element_list[num_packed++] = element;
    }
    element_list.resize(num_packed);
    element_list.push_back(pivot_);
  }
}

inline void QuotientGraph::UnflagPivotStructure() {
  for (const Int& i : elements_[pivot_]) {
    pivot_mask_[i] = 0;
  }
#ifdef QUOTIENT_DEBUG
  for (std::size_t index = 0; index < pivot_mask_.size(); ++index) {
    if (pivot_mask_[index]) {
      std::cerr << "pivot_mask[" << index << "] was " << pivot_mask_[index]
                << " after UnflagPivotStructure." << std::endl;
    }
  }
#endif
}

inline void QuotientGraph::FlagPivotElementList() {
  for (const Int& element : element_lists_[pivot_]) {
    pivot_mask_[element] = -1;
  }
}

inline void QuotientGraph::UnflagPivotElementList() {
  for (const Int& element : element_lists_[pivot_]) {
    pivot_mask_[element] = 0;
  }
#ifdef QUOTIENT_DEBUG
  // At this point, only the principal members of the structure of the
  // current pivot should be flagged (with +1's).
  for (std::size_t index = 0; index < pivot_mask_.size(); ++index) {
    if (pivot_mask_[index] < 0) {
      std::cerr << "pivot_mask[" << index << "] was " << pivot_mask_[index]
                << " after UnflagPivotElementList." << std::endl;
    }
  }
#endif
}

inline void QuotientGraph::AggressiveAbsorption(
    const std::vector<Int>& aggressive_absorption_elements) {
  if (aggressive_absorption_elements.empty()) {
    return;
  }

  for (const Int& e : aggressive_absorption_elements) {
    // Aggressive element absorption:
    //   E_p := (E_p \ E_e) \cup {e}.
    // Instead of explicitly forming this update, we will disable the
    // elements from E_e in the pivot_mask and implicitly keep track of
    // the fact that {e} should be inserted.
    if (control_.store_aggressive_absorptions) {
      aggressive_absorptions_.emplace_back(pivot_, e);
    }
    for (const Int& element : element_lists_[e]) {
      pivot_mask_[element] = 0;
    }
  }

  // Form E_p.
  Int num_packed = 0;
  std::vector<Int>& element_list = element_lists_[pivot_];
  for (const Int& element : element_list) {
    if (pivot_mask_[element]) {
      element_list[num_packed++] = element;
    }
  }
  element_list.resize(num_packed);
  for (const Int& element : aggressive_absorption_elements) {
    if (pivot_mask_[element]) {
      continue;
    }
    element_list.push_back(element);
  }
}

inline Int QuotientGraph::ExactEmptyExternalDegree(Int i) const {
  // Add the cardinality of A_i \ supernode(i).
  Int degree = 0;
  for (const Int& index : adjacency_lists_[i]) {
    degree += supernode_sizes_[index];
  }

#ifdef QUOTIENT_DEBUG
  // We should only have one member of the element list, 'pivot'.
  if (element_lists_[i].size() != 0) {
    std::cerr << "The element list was assumed empty." << std::endl;
  }
#endif

  return degree;
}

inline Int QuotientGraph::ExactSingleExternalDegree(Int i) const {
  // Add the cardinality of A_i \ supernode(i).
  Int degree = 0;
  for (const Int& index : adjacency_lists_[i]) {
    degree += supernode_sizes_[index];
  }

#ifdef QUOTIENT_DEBUG
  // We should only have one member of the element list, 'pivot'.
  if (element_lists_[i].size() != 1) {
    std::cerr << "There was more than one member in the element list."
              << std::endl;
  }
  if (element_lists_[i][0] != pivot_) {
    std::cerr << "The element list should have only contained the pivot."
              << std::endl;
  }
#endif
  // Add |L_p \ supernode(i)|.
  degree += element_sizes_[pivot_] - supernode_sizes_[i];

  return degree;
}

inline Int QuotientGraph::ExactDoubleExternalDegree(Int i) const {
  // Add the cardinality of A_i \ supernode(i).
  Int degree = 0;
  for (const Int& index : adjacency_lists_[i]) {
    degree += supernode_sizes_[index];
  }

#ifdef QUOTIENT_DEBUG
  // We should only have one member of the element list, 'pivot'.
  if (element_lists_[i].size() != 2) {
    std::cerr << "There was more than one member in the element list."
              << std::endl;
  }
  if (element_lists_[i][0] != pivot_ && element_lists_[i][1] != pivot_) {
    std::cerr << "The element list should have contained the pivot."
              << std::endl;
  }
#endif

  // Add |L_p \ supernode(i)|.
  degree += element_sizes_[pivot_] - supernode_sizes_[i];

  // Add |L_e \ L_p|.
  const Int element = element_lists_[i][0] == pivot_ ? element_lists_[i][1] :
      element_lists_[i][0];
  degree += external_element_sizes_[element];

  return degree;
}

inline Int QuotientGraph::ExactGenericExternalDegree(Int i) const {
  // Add the cardinality of A_i \ supernode(i).
  Int degree = 0;
  for (const Int& index : adjacency_lists_[i]) {
    degree += supernode_sizes_[index];
  }

  // Add on the number of unique entries in the structures of the element lists
  // that are outside supernode(i).
  for (const Int& element : element_lists_[i]) {
    for (const Int& j : elements_[element]) {
      if (exact_degree_mask_[j] || i == j || supernode_sizes_[j] < 0) {
        continue;
      }
      degree += supernode_sizes_[j];
      exact_degree_mask_[j] = 1;
    }
  }

  // Clear the mask.
  // NOTE: We allow (*mask)[i] to be set to zero to avoid branching.
  for (const Int& element : element_lists_[i]) {
    for (const Int& j : elements_[element]) {
      exact_degree_mask_[j] = 0;
    }
  }

  return degree;
}

inline Int QuotientGraph::ExactExternalDegree(Int i) const {
  const Int num_elements = element_lists_[i].size();
  if (num_elements == 0) {
    return ExactEmptyExternalDegree(i);
  } else if (num_elements == 1) {
    return ExactSingleExternalDegree(i);
  } else if (num_elements == 2) {
    return ExactDoubleExternalDegree(i);
  } else {
    return ExactGenericExternalDegree(i);
  }
}

inline Int QuotientGraph::AmestoyExternalDegree(Int i) const {
  const Int num_vertices_left =
      num_original_vertices_ - num_eliminated_vertices_;
  const Int old_degree = degree_lists_.degrees[i];

  // Note that this usage of 'external' refers to |L_p \ supernode(i)| and not
  // |L_e \ L_p|, as is the case for 'external_element_sizes'.
  Int external_pivot_structure_size =
      element_sizes_[pivot_] - supernode_sizes_[i];
#ifdef QUOTIENT_DEBUG
  if (external_pivot_structure_size < 0) {
    std::cerr << "Encountered a negative external_pivot_structure_size"
              << std::endl;
  }
#endif

  const Int bound0 = num_vertices_left - supernode_sizes_[i];
  const Int bound1 = old_degree + external_pivot_structure_size;

  // bound_2 = |A_i \ supernode(i)| + |L_p \ supernode(i)| +
  //           \sum_{e in E_i \ {p}} |L_e \ L_p|.
  Int bound2 = 0;
  for (const Int& index : adjacency_lists_[i]) {
#ifdef QUOTIENT_DEBUG
    if (supernode_sizes_[index] < 0) {
      std::cerr << "Encountered an element in an adjacency list." << std::endl;
    }
#endif
    bound2 += supernode_sizes_[index];
  }
  bound2 += external_pivot_structure_size;
  for (const Int& element : element_lists_[i]) {
    if (element == pivot_) {
      continue;
    }
    if (external_element_sizes_[element] >= 0) {
      bound2 += external_element_sizes_[element];
    } else {
#ifdef QUOTIENT_DEBUG
      if (element_sizes_[element] < 0) {
        std::cerr << "Encountered a negative element_size." << std::endl;
      }
#endif
      bound2 += element_sizes_[element];
    }
  }

  const Int degree = std::min(bound0, std::min(bound1, bound2));
  return degree;
}

inline Int QuotientGraph::AshcraftExternalDegree(Int i) const {
  if (element_lists_[i].size() == 2) {
    return ExactDoubleExternalDegree(i);
  }
  return GilbertExternalDegree(i);
}

inline Int QuotientGraph::GilbertExternalDegree(Int i) const {
  Int degree = 0;
  for (const Int& index : adjacency_lists_[i]) {
    degree += supernode_sizes_[index];
  }

  const Int supernode_size = supernode_sizes_[i];
  for (const Int& element : element_lists_[i]) {
    degree += element_sizes_[element] - supernode_size;
  }
  const Int num_vertices_left =
      num_original_vertices_ - num_eliminated_vertices_;
  return std::min(degree, num_vertices_left - supernode_sizes_[i]);
}

inline Int QuotientGraph::ExternalDegree(Int i) const {
  Int degree = -1;
  switch(control_.degree_type) {
    case kExactExternalDegree: {
      degree = ExactExternalDegree(i);
      break;
    }
    case kAmestoyExternalDegree: {
      degree = AmestoyExternalDegree(i);
      break;
    }
    case kAshcraftExternalDegree: {
      degree = AshcraftExternalDegree(i);
      break;
    }
    case kGilbertExternalDegree: {
      degree = GilbertExternalDegree(i);
      break;
    }
  }
  return degree;
}

inline void QuotientGraph::ComputeExternalDegrees(
    std::vector<Int>* external_degrees) {
  const std::size_t supernodal_struct_size = elements_[pivot_].size();
  external_degrees->resize(supernodal_struct_size);
  // TODO(Jack Poulson): Consider how to parallelize using different masks
  // for each thread.
  for (std::size_t index = 0; index < supernodal_struct_size; ++index) {
    const Int i = elements_[pivot_][index];
#ifdef QUOTIENT_DEBUG
    std::vector<Int> degrees(4);
    ExternalDegreeType degree_type = control_.degree_type;
    for (Int type = 0; type < 4; ++type) {
      control_.degree_type = static_cast<ExternalDegreeType>(type);
      degrees[type] = ExternalDegree(i);
    }
    control_.degree_type = degree_type;
    if (degrees[1] < degrees[0] ||
        degrees[2] < degrees[1] ||
        degrees[3] < degrees[2]) {
      std::cerr << "Degrees: " << degrees[0] << ", " << degrees[1] << ", "
                << degrees[2] << ", " << degrees[3] << std::endl;
    }
    (*external_degrees)[index] = degrees[degree_type];
#else
    (*external_degrees)[index] = ExternalDegree(i);
#endif
  }
}

inline void QuotientGraph::UpdateExternalDegrees(
    const std::vector<Int>& external_degrees) {
  const Int supernodal_struct_size = elements_[pivot_].size();
  for (Int index = 0; index < supernodal_struct_size; ++index) {
    const Int i = elements_[pivot_][index];
    const Int degree = external_degrees[index];
    degree_lists_.UpdateDegree(i, degree);
  }
}

inline std::size_t QuotientGraph::VariableHash(
    Int i, VariableHashType hash_type) const {
  if (hash_type == kAshcraftVariableHash) {
    return AshcraftVariableHash(i);
  } else {
    return BasicVariableHash(i);
  }
}

inline std::size_t QuotientGraph::AshcraftVariableHash(Int i) const {
  std::size_t result = BasicVariableHash(i);
  result %= num_original_vertices_;
  return result;
}

inline std::size_t QuotientGraph::BasicVariableHash(Int i) const {
  std::size_t result = 0;
  for (const Int& j : adjacency_lists_[i]) {
    if (supernode_sizes_[j] > 0) {
      result += j;
    }
  }
  for (const Int& j : element_lists_[i]) {
    result += j;
  }
  return result;
}

inline void QuotientGraph::ComputeVariableHashes(
    std::vector<std::size_t>* bucket_keys) {
  const std::size_t supernodal_struct_size = elements_[pivot_].size();
  bucket_keys->resize(supernodal_struct_size);
  #pragma omp parallel for schedule(dynamic)
  for (std::size_t i_index = 0; i_index < supernodal_struct_size; ++i_index) {
    const Int i = elements_[pivot_][i_index];
    (*bucket_keys)[i_index] = VariableHash(i, control_.hash_type);
  }
}

inline bool QuotientGraph::StructuralSupervariablesAreQuotientIndistinguishable(
    Int i, Int j) const {
  if (i == j) {
    return true;
  }
  if (i > j) {
    std::swap(i, j);
  }

  // Check if E_i = E_j.
  if (element_lists_[i].size() != element_lists_[j].size()) {
    return false;
  }
  for (std::size_t index = 0; index < element_lists_[i].size(); ++index) {
    if (element_lists_[i][index] != element_lists_[j][index]) {
      return false;
    }
  }

  // Check if A_i = A_j.
  if (adjacency_lists_[i].size() != adjacency_lists_[j].size()) {
    return false;
  }
  for (std::size_t index = 0; index < adjacency_lists_[i].size(); ++index) {
    if (adjacency_lists_[i][index] != adjacency_lists_[j][index]) {
      return false;
    }
  }

  return true;
}

inline void QuotientGraph::MergeVariables(
    const std::vector<std::size_t>& bucket_keys) {
  // Fill a set of buckets for the hashes of the supernodes adjacent to
  // the current pivot.
  const Int supernodal_struct_size = bucket_keys.size();
  std::vector<Int> bucket_indices;
  for (Int i_index = 0; i_index < supernodal_struct_size; ++i_index) {
    const std::size_t bucket_key = bucket_keys[i_index];
    const Int bucket_index = bucket_key % num_original_vertices_;
    buckets_[bucket_index].push_back(i_index);
    if (buckets_[bucket_index].size() == 2) {
      bucket_indices.push_back(bucket_index);
    }
  }

  std::vector<VariableMergeInfo> variable_merges;
  #pragma omp parallel for schedule(dynamic)
  for (std::size_t index = 0; index < bucket_indices.size(); ++index) {
    const Int bucket_index = bucket_indices[index];
    const std::vector<Int>& bucket = buckets_[bucket_index];
    const Int bucket_size = bucket.size();
    Int num_merges = 0;
    std::vector<bool> merged_supernode(bucket_size, false);
    for (Int i_bucket  = 0; i_bucket < bucket_size; ++i_bucket) {
      if (merged_supernode[i_bucket]) {
        continue;
      }
      const Int i_index = bucket[i_bucket];
      const Int i = elements_[pivot_][i_index];
      for (Int j_bucket = i_bucket + 1; j_bucket < bucket_size; ++j_bucket) {
        if (merged_supernode[j_bucket]) {
          continue;
        }
        const Int j_index = bucket[j_bucket];
        const Int j = elements_[pivot_][j_index];
        if (StructuralSupervariablesAreQuotientIndistinguishable(i, j)) {
          #pragma omp critical
          variable_merges.emplace_back(i, j, supernode_sizes_[j]);

          // Merge [i] -> [j] (and i becomes the principal member).
          next_index_[tail_index_[i]] = j;
          tail_index_[i] = tail_index_[j];
          supernode_sizes_[i] += supernode_sizes_[j];
          supernode_sizes_[j] = 0;

          SwapClearVector(&adjacency_lists_[j]);
          SwapClearVector(&element_lists_[j]);

          ++num_merges;
          merged_supernode[j_bucket] = true;
        }
      }
    }
    num_hash_collisions_ += bucket_size - (num_merges + 1);
  }

  // Update the external degrees in a batched manner.
  for (const VariableMergeInfo& merge : variable_merges) {
    const Int old_degree = degree_lists_.degrees[merge.primary_index];
    const Int degree = old_degree - merge.absorbed_size;
    degree_lists_.UpdateDegree(merge.primary_index, degree);
    degree_lists_.RemoveDegree(merge.absorbed_index);
  }

  // If requested, keep track of the variable merges.
  if (control_.store_variable_merges) {
    for (const VariableMergeInfo& merge : variable_merges) {
      variable_merges_.emplace_back(merge.primary_index, merge.absorbed_index);
    }
  }

  // Clear the buckets.
  for (Int i_index = 0; i_index < supernodal_struct_size; ++i_index) {
    const std::size_t bucket_key = bucket_keys[i_index];
    const Int bucket_index = bucket_key % num_original_vertices_;
    buckets_[bucket_index].clear();
  }
}

inline void QuotientGraph::ConvertPivotIntoElement() {
  const Int supernode_size = supernode_sizes_[pivot_];
  const std::vector<Int>& pivot_element_list = control_.aggressive_absorption ?
      original_pivot_element_list_ : elements_[pivot_];

  // Since this supervariable is being eliminated, it needs to be implicitly
  // removed from all elements containing it. This is accomplished through
  // decrementing their element sizes and marking the supernode as eliminated
  // by making its sign negative.
  for (const Int& element : pivot_element_list) {
    element_sizes_[element] -= supernode_size;
  }
  supernode_sizes_[pivot_] *= -1;

  degree_lists_.RemoveDegree(pivot_);
  SwapClearVector(&adjacency_lists_[pivot_]);
  SwapClearVector(&element_lists_[pivot_]);

  elimination_order_.push_back(pivot_);
  num_eliminated_vertices_ += supernode_size;
}

inline void QuotientGraph::RecomputeExternalElementSizes(
    std::vector<Int>* aggressive_absorption_elements) {
  // Follow the advice at the beginning of Section 5 of [AMD-96] and absorb
  // any element e that satisfies |L_e \ L_p| = 0.
  aggressive_absorption_elements->clear();

  for (const Int& i : elements_[pivot_]) {
    const Int supernode_i_size = supernode_sizes_[i];
#ifdef QUOTIENT_DEBUG
    if (supernode_i_size <= 0) {
      std::cerr << "supernode_" << i << "_size=" << supernode_i_size
                << " in ExternalElementSizes" << std::endl;
    }
#endif
    for (const Int& element : element_lists_[i]) {
      if (element == pivot_) {
        continue;
      }
      Int& external_element_size = external_element_sizes_[element];
      if (external_element_size < 0) {
        external_element_size = element_sizes_[element];
      }

      external_element_size -= supernode_i_size;
      if (control_.aggressive_absorption && external_element_size == 0) {
        aggressive_absorption_elements->push_back(element);
      }
    }
  }
}

inline void QuotientGraph::ResetExternalElementSizes() {
  for (const Int& i : elements_[pivot_]) {
    for (const Int& element : element_lists_[i]) {
      external_element_sizes_[element] = -1;
    }
  }

#ifdef QUOTIENT_DEBUG
  for (std::size_t index = 0; index < external_element_sizes_.size(); ++index) {
    if (external_element_sizes_[index] != -1) {
      std::cerr << "external_element_sizes[" << index << "] was "
                << external_element_sizes_[index] << " after reset."
                << std::endl;
    }
  }
#endif
}

} // namespace quotient

#endif // ifndef QUOTIENT_QUOTIENT_GRAPH_IMPL_H_
