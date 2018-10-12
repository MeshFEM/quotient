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
  supernodes.resize(num_original_vertices);
  for (Int source = 0; source < num_original_vertices; ++source) {
    supernodes[source].push_back(source);
  }

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
    if (supernodes[i].empty()) {
      continue;
    }
    std::cout << "Supernode " << i << "\n";
    PrintVector(supernodes[i], "  members");
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

inline std::vector<Int> QuotientGraph::FormSupernodalStructure(Int element)
    const {
  std::vector<Int> supernodal_structure;
  for (const Int& i : structures[element]) {
    if (!supernodes[i].empty()) {
      supernodal_structure.push_back(i);
    }
  }
  return supernodal_structure;
}

inline std::vector<Int> QuotientGraph::FormSupernodalAdjacencyList(Int i)
    const {
  std::vector<Int> supernodal_adjacency_list;
  for (const Int& j : adjacency_lists[i]) {
    if (!supernodes[j].empty()) {
      supernodal_adjacency_list.push_back(j);
    }
  }
  return supernodal_adjacency_list;
}

inline std::size_t QuotientGraph::AshcraftVariableHash(Int i) const {
  std::size_t result = 0;
  for (const Int& index : adjacency_lists[i]) {
    if (!supernodes[index].empty()) {
      result += index;
      result %= num_original_vertices - 1;
    }
  }
  for (const Int& index : element_lists[i]) {
    result += index;
    result %= num_original_vertices - 1;
  }
  return result + 1;
};

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
    const std::vector<Int>& supernodal_structure,
    bool aggressive_absorption,
    std::vector<Int>* external_structure_sizes,
    std::vector<Int>* aggressive_absorption_elements) const {
  // Follow the advice at the beginning of Section 5 of [AMD-96] and absorb
  // any element e that satisfies |L_e \ L_p| = 0.
  aggressive_absorption_elements->clear();

  for (const Int& i : supernodal_structure) {
    for (const Int& element : element_lists[i]) {
      if ((*external_structure_sizes)[element] < 0) {
        (*external_structure_sizes)[element] = structures[element].size();
      }
      Int& external_structure_size = (*external_structure_sizes)[element];
      external_structure_size -= supernodes[i].size();
      if (aggressive_absorption && external_structure_size == 0) {
        aggressive_absorption_elements->push_back(element);
      }
#ifdef QUOTIENT_DEBUG
      for (const Int& j : supernodes[i]) {
        auto iter = std::lower_bound(
            structures[element].begin(), structures[element].end(), j);
        if (iter == structures[element].end() || *iter != j) {
          std::cerr << "structures[" << element << "] did not contain a full "
                       "copy of supernode " << i << std::endl;
        }
      }
#endif
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
