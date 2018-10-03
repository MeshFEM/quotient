/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_QUOTIENT_GRAPH_H_
#define QUOTIENT_QUOTIENT_GRAPH_H_

#include <iostream>
#include <set>
#include <vector>

#include "quotient/config.hpp"
#include "quotient/coordinate_graph.hpp"
#include "quotient/random_access_heap.hpp"

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
struct QuotientGraph {
  // The number of vertices in the original graph.
  Int num_original_vertices;

  // The number of vertices that have been eliminated from the original graph.
  Int num_eliminated_vertices;

  // A list of length 'num_original_vertices' of supernodal sets. If index
  // 'i' is a principal vertex of a supernode, then 'supernodes[i]'
  // is a nontrivial, sorted list of integers that composes supernode(i).
  std::vector<std::vector<Int>> supernodes;

  // A list of length 'num_original_vertices' of element nonzero structures.
  // The 'element' index of the list, 'structures[element]', will be created
  // when supernode 'element' is converted from a variable to an element.
  //
  // The structure lists connect to each individual member of any supernode.
  std::vector<std::vector<Int>> structures;

  // A list of length 'num_original_vertices' of the (unmodified) variable
  // adjacencies of each principal variable. For example, if index 'i' is a
  // principal variable, then 'adjacency_lists[i]' contains the set of neighbor
  // variables for variable i that are not redundant with respect to edges
  // implied by 'structures'.
  //
  // The adjacency lists connect to each individual member of any supernode.
  std::vector<std::vector<Int>> adjacency_lists;

  // A list of length 'num_original_vertices' of the elements adjacent to
  // each principal variable. For example, if index 'i' is a principal
  // variable, then 'element_lists[i]' contains the list of elements adjacent
  // to supervariable 'i'.
  //
  // The element lists only contain the principal member of any supernode.
  std::vector<std::vector<Int>> element_lists;

#ifdef QUOTIENT_DEBUG
  // The list of all (principal) variables.
  std::set<Int> variables;

  // The list of all (principal) elements (principal members of eliminated
  // variables).
  std::set<Int> elements;
#endif

  // Initializes the quotient graph from a symmetric graph.
  QuotientGraph(const CoordinateGraph& graph);

  // Uses 'supernodes' to filter 'structures[element]' into just its entries
  // which are principal members of a supernode.
  std::vector<Int> FormSupernodalStructure(Int element) const;

  // Uses 'supernodes' to filter 'adjaency_lists[i]' into just its entries
  // which are principal members of a supernode.
  std::vector<Int> FormSupernodalAdjacencyList(Int i) const;
};

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
      if (edge.target != source) {
        adjacency_list.push_back(edge.target);
      }
    }
  }

  // Trivially initialize the element lists.
  element_lists.resize(num_original_vertices);

#ifdef QUOTIENT_DEBUG
  // Initialize the (principal) variable and element lists.
  for (Int source = 0; source < num_original_vertices; ++source) {
    variables.insert(source);
  }
#endif
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

// Pretty-prints an std::vector<T>.
template<typename T>
void PrintVector(const std::vector<T>& vec, const std::string& msg) {
  std::cout << msg << ": ";
  for (std::size_t i = 0; i < vec.size(); ++i) {
    std::cout << vec[i] << " ";
  }
  std::cout << "\n";
}

// Pretty-prints a QuotientGraph.
void PrintGraph(
    const QuotientGraph& graph,
    const RandomAccessHeap<Int>& external_degree_heap) {
  for (Int i = 0; i < graph.num_original_vertices; ++i) {
    if (graph.supernodes[i].empty()) {
      continue;
    }
    std::cout << "Supernode " << i << "\n";
    PrintVector(graph.supernodes[i], "  members");
    std::cout << "  external_degree: " << external_degree_heap.Value(i) << "\n";
    if (external_degree_heap.ValidValue(i)) {
      PrintVector(graph.adjacency_lists[i], "  adjacency_list");
      PrintVector(graph.element_lists[i], "  element_list");
    } else {
      PrintVector(graph.structures[i], "  structure");
    }
    std::cout << "\n";
  }
}

} // namespace quotient

#endif // ifndef QUOTIENT_QUOTIENT_GRAPH_H_
