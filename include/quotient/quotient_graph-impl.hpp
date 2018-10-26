/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_QUOTIENT_GRAPH_IMPL_H_
#define QUOTIENT_QUOTIENT_GRAPH_IMPL_H_

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

inline QuotientGraph::QuotientGraph()
: num_original_vertices(0), num_eliminated_vertices(0) { }

inline QuotientGraph::QuotientGraph(
    const CoordinateGraph& graph,
    const MinimumDegreeControl& control)
: num_original_vertices(graph.NumSources()), num_eliminated_vertices(0),
  control_(control) {
  // Initialize the supernodes as simple.
  supernode_sizes.resize(num_original_vertices, 1);
  next_index.resize(num_original_vertices, -1);
  tail_index.resize(num_original_vertices);
  std::iota(tail_index.begin(), tail_index.end(), 0);

  // Initialize the adjacency lists from the original graph.
  adjacency_lists.resize(num_original_vertices);
  const std::vector<GraphEdge>& edges = graph.Edges();
  for (Int source = 0; source < num_original_vertices; ++source) {
    const Int source_edge_offset = graph.SourceEdgeOffset(source);
    const Int next_source_edge_offset = graph.SourceEdgeOffset(source + 1);
    const Int num_connections = next_source_edge_offset - source_edge_offset;
    std::vector<Int>& adjacency_list = adjacency_lists[source];

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
  degree_lists.degrees.resize(num_original_vertices);
  degree_lists.degree_heads.resize(num_original_vertices - 1, -1);
  degree_lists.next_degree_member.resize(num_original_vertices, -1);
  degree_lists.last_degree_member.resize(num_original_vertices, -1);
  for (Int source = 0; source < num_original_vertices; ++source) {
    const Int degree = adjacency_lists[source].size();
    degree_lists.AddDegree(source, degree);
  }

  // Trivially initialize the element lists.
  element_lists.resize(num_original_vertices);

  // Trivially initialize the lower-triangular nonzero structures.
  if (control_.store_structures) {
    structures.resize(num_original_vertices);
  }
  elements.resize(num_original_vertices);
  element_sizes.resize(num_original_vertices, 0);

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
    external_element_sizes_.resize(num_original_vertices, -1);
  }

  pivot_mask_.resize(num_original_vertices, 0);

#ifdef QUOTIENT_DEBUG
  exact_degree_mask_.resize(num_original_vertices, 0);
#else
  if (control_.degree_type == kExactExternalDegree) {
    exact_degree_mask_.resize(num_original_vertices, 0);
  }
#endif
}

inline Int QuotientGraph::GetNextPivot() {
  pivot_ = degree_lists.FindMinimalIndex(control_.force_minimal_pivot_indices);
  return pivot_;
}

inline bool QuotientGraph::UsingExternalElementSizes() const {
  return using_external_element_sizes_;
}

inline const std::vector<std::pair<Int, Int>>&
QuotientGraph::VariableMerges() const {
  return variable_merges_;
}

inline void QuotientGraph::Print() const {
  for (Int i = 0; i < num_original_vertices; ++i) {
    if (supernode_sizes[i] == 0) {
      continue;
    }
    std::cout << "Supernode " << i << "\n";
    const std::vector<Int> supernode = FormSupernode(i);
    PrintVector(supernode, "  members");
    if (supernode_sizes[i] < 0) {
      PrintVector(elements[i], "  elements");
    } else {
      PrintVector(adjacency_lists[i], "  adjacency_list");
      PrintVector(element_lists[i], "  element_list");
    }
    std::cout << "\n";
  }
}

inline std::vector<Int> QuotientGraph::FormSupernode(Int i) const {
  Int index = i;
  Int num_members = 1;
  while (index != tail_index[i]) {
    index = next_index[index];
    ++num_members;
  }

  std::vector<Int> supernode;
  supernode.reserve(num_members);
  index = i;
  supernode.push_back(index);
  while (index != tail_index[i]) {
    index = next_index[index];
    supernode.push_back(index);
  }

  return supernode;
}

inline Int QuotientGraph::ComputePivotStructure() {
#ifdef QUOTIENT_DEBUG
  for (Int index = 0; index < num_original_vertices; ++index) {
    if (pivot_mask_[index]) {
      std::cerr << "Pivot mask had a nonzero on entry to ComputePivotStructure."
                << std::endl;
    }
  }
  if (!elements[pivot_].empty()) {
    std::cerr << "Chose a pivot more than once." << std::endl;
  }
#endif
  // TODO(Jack Poulson): Reserve space for graph->elements[pivot].
  for (const Int& index : adjacency_lists[pivot_]) {
    const Int supernode_size = supernode_sizes[index];
#ifdef QUOTIENT_DEBUG
    if (supernode_size < 0) {
      std::cerr << "Encountered an element in the adjacency list." << std::endl;
    }
#endif
    if (!supernode_size) {
      continue;
    }
    pivot_mask_[index] = 1;
    elements[pivot_].push_back(index);
    element_sizes[pivot_] += supernode_size;

    if (control_.store_structures) {
      // Fill in the supernode.
      Int k = index;
      structures[pivot_].push_back(k);
      while (k != tail_index[index]) {
        k = next_index[k];
        structures[pivot_].push_back(k);
      }
    }
  }

  Int num_stale_element_members = 0;
  for (const Int& element : element_lists[pivot_]) {
    for (const Int& index : elements[element]) {
      const Int supernode_size = supernode_sizes[index];
      if (supernode_size < 0) {
        ++num_stale_element_members;
        // TODO(Jack Poulson): Quickly remove elements[element][index]. Though
        // this branch appears to be extremely rare.
      }
      if (index == pivot_ || supernode_size <= 0 || pivot_mask_[index]) {
        continue;
      }
      pivot_mask_[index] = 1;
      elements[pivot_].push_back(index);
      element_sizes[pivot_] += supernode_size;

      if (control_.store_structures) {
        // Fill in the supernode.
        Int k = index;
        structures[pivot_].push_back(k);
        while (k != tail_index[index]) {
          k = next_index[k];
          structures[pivot_].push_back(k);
        }
      }
    }
  }

  // If aggressive absorption is to be activated, we need to save a copy of
  // the pivot element list so that element_sizes can be later decremented after
  // converting the pivot into an element.
  if (control_.aggressive_absorption) {
    original_pivot_element_list_ = element_lists[pivot_];
  }

#ifdef QUOTIENT_DEBUG
  for (const Int& element : element_lists[pivot_]) {
    Int element_size = 0;
    for (const Int& index : elements[element]) {
      const Int supernode_size = supernode_sizes[index];
      if (supernode_size <= 0) {
        continue;
      }
      element_size += supernode_size;
    }
    const Int cached_element_size = element_sizes[element];
    if (element_size != cached_element_size) {
      std::cerr << "Element " << element << " had a size of "
                << element_size << " but the cached size was "
                << cached_element_size << std::endl;
    }
  }
#endif

  return num_stale_element_members;
}

inline void QuotientGraph::UpdateAdjacencyListsAfterSelectingPivot() {
  for (const Int& i : elements[pivot_]) {
    // Remove redundant adjacency entries:
    //   A_i := (A_i \ L_p) \ supernode(p).
    Int packed_size = 0;
    for (Int index : adjacency_lists[i]) {
      if (index == pivot_ || pivot_mask_[index] || !supernode_sizes[index]) {
        continue;
      }
      adjacency_lists[i][packed_size++] = index;
    }
    adjacency_lists[i].resize(packed_size);
  }
}

inline void QuotientGraph::UpdateElementListsAfterSelectingPivot() {
  for (const Int& i : elements[pivot_]) {
    // Element absorption:
    //   E_i := (E_i \ E_p) \cup {p}.
    std::vector<Int>& element_list = element_lists[i];
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
  for (const Int& i : elements[pivot_]) {
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
  for (const Int& element : element_lists[pivot_]) {
    pivot_mask_[element] = -1;
  }
}

inline void QuotientGraph::UnflagPivotElementList() {
  for (const Int& element : element_lists[pivot_]) {
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
      aggressive_absorptions.emplace_back(pivot_, e);
    }
    for (const Int& element : element_lists[e]) {
      pivot_mask_[element] = 0;
    }
  }

  // Form E_p.
  Int num_packed = 0;
  std::vector<Int>& element_list = element_lists[pivot_];
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
  for (const Int& index : adjacency_lists[i]) {
    degree += supernode_sizes[index];
  }

#ifdef QUOTIENT_DEBUG
  // We should only have one member of the element list, 'pivot'.
  if (element_lists[i].size() != 0) {
    std::cerr << "The element list was assumed empty." << std::endl;
  }
#endif

  return degree;
}

inline Int QuotientGraph::ExactSingleExternalDegree(Int i) const {
  // Add the cardinality of A_i \ supernode(i).
  Int degree = 0;
  for (const Int& index : adjacency_lists[i]) {
    degree += supernode_sizes[index];
  }

#ifdef QUOTIENT_DEBUG
  // We should only have one member of the element list, 'pivot'.
  if (element_lists[i].size() != 1) {
    std::cerr << "There was more than one member in the element list."
              << std::endl;
  }
  if (element_lists[i][0] != pivot_) {
    std::cerr << "The element list should have only contained the pivot."
              << std::endl;
  }
#endif
  // Add |L_p \ supernode(i)|.
  degree += element_sizes[pivot_] - supernode_sizes[i];

  return degree;
}

inline Int QuotientGraph::ExactDoubleExternalDegree(Int i) const {
  // Add the cardinality of A_i \ supernode(i).
  Int degree = 0;
  for (const Int& index : adjacency_lists[i]) {
    degree += supernode_sizes[index];
  }

#ifdef QUOTIENT_DEBUG
  // We should only have one member of the element list, 'pivot'.
  if (element_lists[i].size() != 2) {
    std::cerr << "There was more than one member in the element list."
              << std::endl;
  }
  if (element_lists[i][0] != pivot_ && element_lists[i][1] != pivot_) {
    std::cerr << "The element list should have contained the pivot."
              << std::endl;
  }
#endif

  // Add |L_p \ supernode(i)|.
  degree += element_sizes[pivot_] - supernode_sizes[i];

  // Add |L_e \ L_p|.
  const Int element = element_lists[i][0] == pivot_ ? element_lists[i][1] :
      element_lists[i][0];
  degree += external_element_sizes_[element];

  return degree;
}

inline Int QuotientGraph::ExactGenericExternalDegree(Int i) const {
  // Add the cardinality of A_i \ supernode(i).
  Int degree = 0;
  for (const Int& index : adjacency_lists[i]) {
    degree += supernode_sizes[index];
  }

  // Add on the number of unique entries in the structures of the element lists
  // that are outside supernode(i).
  for (const Int& element : element_lists[i]) {
    for (const Int& j : elements[element]) {
      if (exact_degree_mask_[j] || i == j || supernode_sizes[j] < 0) {
        continue;
      }
      degree += supernode_sizes[j];
      exact_degree_mask_[j] = 1;
    }
  }

  // Clear the mask.
  // NOTE: We allow (*mask)[i] to be set to zero to avoid branching.
  for (const Int& element : element_lists[i]) {
    for (const Int& j : elements[element]) {
      exact_degree_mask_[j] = 0;
    }
  }

  return degree;
}

inline Int QuotientGraph::ExactExternalDegree(Int i) const {
  const Int num_elements = element_lists[i].size();
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
  const Int num_vertices_left = num_original_vertices - num_eliminated_vertices;
  const Int old_degree = degree_lists.degrees[i];

  // Note that this usage of 'external' refers to |L_p \ supernode(i)| and not
  // |L_e \ L_p|, as is the case for 'external_element_sizes'.
  Int external_pivot_structure_size =
      element_sizes[pivot_] - supernode_sizes[i];
#ifdef QUOTIENT_DEBUG
  if (external_pivot_structure_size < 0) {
    std::cerr << "Encountered a negative external_pivot_structure_size"
              << std::endl;
  }
#endif

  const Int bound0 = num_vertices_left - supernode_sizes[i];
  const Int bound1 = old_degree + external_pivot_structure_size;

  // bound_2 = |A_i \ supernode(i)| + |L_p \ supernode(i)| +
  //           \sum_{e in E_i \ {p}} |L_e \ L_p|.
  Int bound2 = 0;
  for (const Int& index : adjacency_lists[i]) {
#ifdef QUOTIENT_DEBUG
    if (supernode_sizes[index] < 0) {
      std::cerr << "Encountered an element in an adjacency list." << std::endl;
    }
#endif
    bound2 += supernode_sizes[index];
  }
  bound2 += external_pivot_structure_size;
  for (const Int& element : element_lists[i]) {
    if (element == pivot_) {
      continue;
    }
    if (external_element_sizes_[element] >= 0) {
      bound2 += external_element_sizes_[element];
    } else {
#ifdef QUOTIENT_DEBUG
      if (element_sizes[element] < 0) {
        std::cerr << "Encountered a negative element_size." << std::endl;
      }
#endif
      bound2 += element_sizes[element];
    }
  }

  const Int degree = std::min(bound0, std::min(bound1, bound2));
  return degree;
}

inline Int QuotientGraph::AshcraftExternalDegree(Int i) const {
  if (element_lists[i].size() == 2) {
    return ExactDoubleExternalDegree(i);
  }
  return GilbertExternalDegree(i);
}

inline Int QuotientGraph::GilbertExternalDegree(Int i) const {
  Int degree = 0;
  for (const Int& index : adjacency_lists[i]) {
    degree += supernode_sizes[index];
  }

  const Int supernode_size = supernode_sizes[i];
  for (const Int& element : element_lists[i]) {
    degree += element_sizes[element] - supernode_size;
  }
  const Int num_vertices_left = num_original_vertices - num_eliminated_vertices;
  return std::min(degree, num_vertices_left - supernode_sizes[i]);
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
  const std::size_t supernodal_struct_size = elements[pivot_].size();
  external_degrees->resize(supernodal_struct_size);
  // TODO(Jack Poulson): Consider how to parallelize using different masks
  // for each thread.
  for (std::size_t index = 0; index < supernodal_struct_size; ++index) {
    const Int i = elements[pivot_][index];
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
  const Int supernodal_struct_size = elements[pivot_].size();
  for (Int index = 0; index < supernodal_struct_size; ++index) {
    const Int i = elements[pivot_][index];
    const Int degree = external_degrees[index];
    degree_lists.UpdateDegree(i, degree);
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
  result %= num_original_vertices - 1;
  return result + 1;
}

inline std::size_t QuotientGraph::BasicVariableHash(Int i) const {
  std::size_t result = 0;
  for (const Int& j : adjacency_lists[i]) {
    if (supernode_sizes[j] > 0) {
      result += j;
    }
  }
  for (const Int& j : element_lists[i]) {
    result += j;
  }
  return result;
}

inline void QuotientGraph::ComputeVariableHashes(
    std::vector<std::size_t>* bucket_keys) {
  const std::size_t supernodal_struct_size = elements[pivot_].size();
  bucket_keys->resize(supernodal_struct_size);
  #pragma omp parallel for schedule(dynamic)
  for (std::size_t i_index = 0; i_index < supernodal_struct_size; ++i_index) {
    const Int i = elements[pivot_][i_index];
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
  if (element_lists[i].size() != element_lists[j].size()) {
    return false;
  }
  for (std::size_t index = 0; index < element_lists[i].size(); ++index) {
    if (element_lists[i][index] != element_lists[j][index]) {
      return false;
    }
  }

  // Check if A_i = A_j.
  if (adjacency_lists[i].size() != adjacency_lists[j].size()) {
    return false;
  }
  for (std::size_t index = 0; index < adjacency_lists[i].size(); ++index) {
    if (adjacency_lists[i][index] != adjacency_lists[j][index]) {
      return false;
    }
  }

  return true;
}

inline void QuotientGraph::MergeVariables(
    const std::vector<std::size_t>& bucket_keys,
    std::vector<std::vector<Int>>* buckets) {
  // Fill a set of buckets for the hashes of the supernodes adjacent to
  // the current pivot.
  const Int supernodal_struct_size = bucket_keys.size();
  std::vector<Int> bucket_indices;
  for (Int i_index = 0; i_index < supernodal_struct_size; ++i_index) {
    const std::size_t bucket_key = bucket_keys[i_index];
    const Int bucket_index = bucket_key % num_original_vertices;
    (*buckets)[bucket_index].push_back(i_index);
    if ((*buckets)[bucket_index].size() == 2) {
      bucket_indices.push_back(bucket_index);
    }
  }

  std::vector<VariableMergeInfo> variable_merges;
  #pragma omp parallel for schedule(dynamic)
  for (std::size_t index = 0; index < bucket_indices.size(); ++index) {
    const Int bucket_index = bucket_indices[index];
    const std::vector<Int>& bucket = (*buckets)[bucket_index];
    const Int bucket_size = bucket.size();
    std::vector<bool> merged_supernode(bucket_size, false);
    for (Int i_bucket  = 0; i_bucket < bucket_size; ++i_bucket) {
      if (merged_supernode[i_bucket]) {
        continue;
      }
      const Int i_index = bucket[i_bucket];
      const Int i = elements[pivot_][i_index];
      for (Int j_bucket = i_bucket + 1; j_bucket < bucket_size; ++j_bucket) {
        if (merged_supernode[j_bucket]) {
          continue;
        }
        const Int j_index = bucket[j_bucket];
        const Int j = elements[pivot_][j_index];
        if (StructuralSupervariablesAreQuotientIndistinguishable(i, j)) {
          #pragma omp critical
          variable_merges.emplace_back(i, j, supernode_sizes[j]);

          // Merge [i] -> [j] (and i becomes the principal member).
          next_index[tail_index[i]] = j;
          tail_index[i] = tail_index[j];
          supernode_sizes[i] += supernode_sizes[j];
          supernode_sizes[j] = 0;

          SwapClearVector(&adjacency_lists[j]);
          SwapClearVector(&element_lists[j]);

          merged_supernode[j_bucket] = true;
        }
      }
    }
  }

  // Update the external degrees in a batched manner.
  for (const VariableMergeInfo& merge : variable_merges) {
    const Int old_degree = degree_lists.degrees[merge.primary_index];
    const Int degree = old_degree - merge.absorbed_size;
    degree_lists.UpdateDegree(merge.primary_index, degree);
    degree_lists.RemoveDegree(merge.absorbed_index);
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
    const Int bucket_index = bucket_key % num_original_vertices;
    (*buckets)[bucket_index].clear();
  }
}

inline void QuotientGraph::ConvertPivotIntoElement() {
  const Int supernode_size = supernode_sizes[pivot_];
  const std::vector<Int>& pivot_element_list = control_.aggressive_absorption ?
      original_pivot_element_list_ : elements[pivot_];

  // Since this supervariable is being eliminated, it needs to be implicitly
  // removed from all elements containing it. This is accomplished through
  // decrementing their element sizes and marking the supernode as eliminated
  // by making its sign negative.
  for (const Int& element : pivot_element_list) {
    element_sizes[element] -= supernode_size;
  }
  supernode_sizes[pivot_] *= -1;

  degree_lists.RemoveDegree(pivot_);
  SwapClearVector(&adjacency_lists[pivot_]);
  SwapClearVector(&element_lists[pivot_]);

  num_eliminated_vertices += supernode_size;
}

inline void QuotientGraph::RecomputeExternalElementSizes(
    std::vector<Int>* aggressive_absorption_elements) {
  // Follow the advice at the beginning of Section 5 of [AMD-96] and absorb
  // any element e that satisfies |L_e \ L_p| = 0.
  aggressive_absorption_elements->clear();

  for (const Int& i : elements[pivot_]) {
    const Int supernode_i_size = supernode_sizes[i];
#ifdef QUOTIENT_DEBUG
    if (supernode_i_size <= 0) {
      std::cerr << "supernode_" << i << "_size=" << supernode_i_size
                << " in ExternalElementSizes" << std::endl;
    }
#endif
    for (const Int& element : element_lists[i]) {
      if (element == pivot_) {
        continue;
      }
      Int& external_element_size = external_element_sizes_[element];
      if (external_element_size < 0) {
        external_element_size = element_sizes[element];
      }

      external_element_size -= supernode_i_size;
      if (control_.aggressive_absorption && external_element_size == 0) {
        aggressive_absorption_elements->push_back(element);
      }
    }
  }
}

inline void QuotientGraph::ResetExternalElementSizes() {
  for (const Int& i : elements[pivot_]) {
    for (const Int& element : element_lists[i]) {
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
