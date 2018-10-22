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

inline QuotientGraph::QuotientGraph(const CoordinateGraph& graph)
: num_original_vertices(graph.NumSources()), num_eliminated_vertices(0) {
  // Initialize the supernodes as simple.
  supernode_sizes.resize(num_original_vertices, 1);
  next_index.resize(num_original_vertices, -1);
  head_index.resize(num_original_vertices);
  std::iota(head_index.begin(), head_index.end(), 0);
  tail_index.resize(num_original_vertices);
  std::iota(tail_index.begin(), tail_index.end(), 0);

  // Trivially initialize the lower-triangular nonzero structures.
  structures.resize(num_original_vertices);

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

  // Trivially initialize the element lists.
  element_lists.resize(num_original_vertices);

  // Initialize the cached binary tree of external degrees.
  std::vector<Int> external_degrees_vec(num_original_vertices);
  for (Int source = 0; source < num_original_vertices; ++source) {
    external_degrees_vec[source] = adjacency_lists[source].size();
  }
  external_degree_heap.Reset(external_degrees_vec);
}

inline void QuotientGraph::Print() const {
  for (Int i = 0; i < num_original_vertices; ++i) {
    if (supernode_sizes[i] <= 0) {
      continue;
    }
    std::cout << "Supernode " << i << "\n";
    const std::vector<Int> supernode = FormSupernode(i);
    PrintVector(supernode, "  members");
    std::cout << "  external_degree: " << external_degree_heap.Value(i) << "\n";
    if (external_degree_heap.ValidValue(i)) {
      PrintVector(adjacency_lists[i], "  adjacency_list");
      PrintVector(element_lists[i], "  element_list");
    } else {
      PrintVector(structures[i], "  structure");
    }
    std::cout << "\n";
  }
}

inline std::vector<Int> QuotientGraph::FormSupernode(Int i) const {
  Int num_members = 0;
  Int index = i;
  while (true) {
    ++num_members;
    if (index == tail_index[i]) {
      break;
    }
    index = next_index[index];
  }

  std::vector<Int> supernode;
  supernode.reserve(num_members);
  index = i;
  while (true) {
    supernode.push_back(index);
    if (index == tail_index[i]) {
      break;
    }
    index = next_index[index];
  }

  return supernode;
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
  std::size_t result = 0;
  for (const Int& j : adjacency_lists[i]) {
    if (supernode_sizes[j] > 0) {
      result += j;
      result %= num_original_vertices - 1;
    }
  }
  for (const Int& j : element_lists[i]) {
    result += j;
    result %= num_original_vertices - 1;
  }
  return result + 1;
}

inline std::size_t QuotientGraph::BasicVariableHash(Int i) const {
  std::size_t result = 0;
  for (std::size_t index = 0; index < adjacency_lists[i].size(); ++index) {
    const Int j = adjacency_lists[i][index];
    result += (index + 1) * (j + 1);
  }
  for (std::size_t index = 0; index < element_lists[i].size(); ++index) {
    const Int j = element_lists[i][index];
    result += (index + 1) * (j + 1);
  }
  return result;
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

inline void QuotientGraph::InitializeExternalStructureSizes(
    std::vector<Int>* external_structure_sizes) const {
  external_structure_sizes->clear();
  external_structure_sizes->resize(num_original_vertices, -1);
}

inline void QuotientGraph::ExternalStructureSizes(
    const std::vector<Int>& supernodal_pivot_structure,
    bool aggressive_absorption,
    std::vector<Int>* external_structure_sizes,
    std::vector<Int>* aggressive_absorption_elements) const {
  // Follow the advice at the beginning of Section 5 of [AMD-96] and absorb
  // any element e that satisfies |L_e \ L_p| = 0.
  aggressive_absorption_elements->clear();

  for (const Int& i : supernodal_pivot_structure) {
    for (const Int& element : element_lists[i]) {
      if ((*external_structure_sizes)[element] < 0) {
        (*external_structure_sizes)[element] = structures[element].size();
      }
      Int& external_structure_size = (*external_structure_sizes)[element];
      external_structure_size -= supernode_sizes[i];
      if (aggressive_absorption && external_structure_size == 0) {
        aggressive_absorption_elements->push_back(element);
      }
    }
  }
}

inline void QuotientGraph::ResetExternalStructureSizes(
    const std::vector<Int>& supernodal_structure,
    std::vector<Int>* external_structure_sizes) const {
  for (const Int& i : supernodal_structure) {
    for (const Int& element : element_lists[i]) {
      (*external_structure_sizes)[element] = -1;
    }
  }
}

} // namespace quotient

#endif // ifndef QUOTIENT_QUOTIENT_GRAPH_IMPL_H_
