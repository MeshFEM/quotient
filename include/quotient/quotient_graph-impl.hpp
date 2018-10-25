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

inline Int DegreeLists::FindMinimalIndex(bool demand_smallest_index) {
  while (degree_lower_bound < Int(degree_heads.size()) &&
      degree_heads[degree_lower_bound] == -1) {
    ++degree_lower_bound;
  }
  if (degree_lower_bound == Int(degree_heads.size())) {
    return -1;
  }

  Int index = degree_heads[degree_lower_bound];
  if (demand_smallest_index) {
    Int minimal_index = index;
    while (next_degree_member[index] != -1) {
      index = next_degree_member[index];
      minimal_index = std::min(minimal_index, index);
    }
    index = minimal_index;
  }
  return index;
}

inline void DegreeLists::RemoveDegree(Int index) {
  const Int degree = degrees[index];
  const Int last = last_degree_member[index];
  const Int next = next_degree_member[index];
  if (last == -1) {
    degree_heads[degree] = next;
  } else {
    // 'index' was not the head, so simply patch the connections.
    next_degree_member[last] = next;
  }
  if (next != -1) {
    last_degree_member[next] = last;
  }
}

inline void DegreeLists::AddDegree(Int index, Int degree) {
  const Int head = degree_heads[degree];
  degree_heads[degree] = index;
  last_degree_member[index] = -1;
  next_degree_member[index] = head;
  if (head != -1) {
    last_degree_member[head] = index;
  }
  degrees[index] = degree;

  // Update the minimal degree.
  degree_lower_bound = std::min(degree_lower_bound, degree);
}

inline void DegreeLists::UpdateDegree(Int index, Int degree) {
  if (degrees[index] == degree) {
    return;
  }
  RemoveDegree(index);
  AddDegree(index, degree);
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
  // TODO(Jack Poulson): Avoid allocating 'structures' unless the
  // 'store_structures' flag of 'MinimumDegreeControl' is true.
  structures.resize(num_original_vertices);
  elements.resize(num_original_vertices);
  element_sizes.resize(num_original_vertices, 0);
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

inline void QuotientGraph::InitializeExternalElementSizes(
    std::vector<Int>* external_element_sizes) const {
  external_element_sizes->clear();
  external_element_sizes->resize(num_original_vertices, -1);
}

inline void QuotientGraph::ExternalElementSizes(
    Int pivot,
    bool aggressive_absorption,
    std::vector<Int>* external_element_sizes,
    std::vector<Int>* aggressive_absorption_elements) const {
  // Follow the advice at the beginning of Section 5 of [AMD-96] and absorb
  // any element e that satisfies |L_e \ L_p| = 0.
  aggressive_absorption_elements->clear();

  for (const Int& i : elements[pivot]) {
    const Int supernode_i_size = supernode_sizes[i];
#ifdef QUOTIENT_DEBUG
    if (supernode_i_size <= 0) {
      std::cerr << "supernode_" << i << "_size=" << supernode_i_size
                << " in ExternalElementSizes" << std::endl;
    }
#endif
    for (const Int& element : element_lists[i]) {
      if (element == pivot) {
        continue;
      }
      Int& external_element_size = (*external_element_sizes)[element];
      if (external_element_size < 0) {
        external_element_size = element_sizes[element];
      }

      external_element_size -= supernode_i_size;
      if (aggressive_absorption && external_element_size == 0) {
        aggressive_absorption_elements->push_back(element);
      }
    }
  }
}

inline void QuotientGraph::ResetExternalElementSizes(
    const std::vector<Int>& original_pivot_element,
    std::vector<Int>* external_element_sizes) const {
  for (const Int& i : original_pivot_element) {
    for (const Int& element : element_lists[i]) {
      (*external_element_sizes)[element] = -1;
    }
  }

#ifdef QUOTIENT_DEBUG
  for (std::size_t index = 0; index < external_element_sizes->size(); ++index) {
    if ((*external_element_sizes)[index] != -1) {
      std::cerr << "external_element_sizes[" << index << "] was "
                << (*external_element_sizes)[index] << " after reset."
                << std::endl;
    }
  }
#endif
}

} // namespace quotient

#endif // ifndef QUOTIENT_QUOTIENT_GRAPH_IMPL_H_
