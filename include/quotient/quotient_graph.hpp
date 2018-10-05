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

// Pretty-prints an std::vector<T>.
// TODO(Jack Poulson): Find a better location for this utility function.
template<typename T>
void PrintVector(const std::vector<T>& vec, const std::string& msg) {
  std::cout << msg << ": ";
  for (std::size_t i = 0; i < vec.size(); ++i) {
    std::cout << vec[i] << " ";
  }
  std::cout << "\n";
}

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

  // A cached binary tree of external degrees that allows for O(lg(n))
  // random modification, O(1) random access, and O(1) extraction of the
  // left-most minimal index.
  RandomAccessHeap<Int> external_degree_heap;

  // An optional list of aggressive element absorption pairs: each pair (e, f)
  // consists of the absorbing element, e, and the absorbed element, f.
  std::vector<std::pair<Int, Int>> aggressive_absorptions;

  // An optional list of supervariable merge pairs: each pair (i, j) consists of
  // the absorbing supervariable, i, and the absorbed supervariable, j.
  std::vector<std::pair<Int, Int>> variable_merges;

  // Initializes the quotient graph from a symmetric graph.
  QuotientGraph(const CoordinateGraph& graph);

  // Pretty-prints the QuotientGraph.
  void Print() const;

  // Uses 'supernodes' to filter 'structures[element]' into just its entries
  // which are principal members of a supernode.
  std::vector<Int> FormSupernodalStructure(Int element) const;

  // Uses 'supernodes' to filter 'adjaency_lists[i]' into just its entries
  // which are principal members of a supernode.
  std::vector<Int> FormSupernodalAdjacencyList(Int i) const;

  // A definition of Ashcraft's hash function (as described in [ADD-96]).
  // Note that only principal members of A_i are incorporated in the hash.
  std::size_t AshcraftVariableHash(Int i) const;

  // Returns true if supernodes 'i' and 'j' are considered indistinguishable
  // with respect to their quotient graph representation.
  //
  // The elimination graph definition (e.g., as given by [ADD-96])
  // involves testing if
  //
  //   Adj_{GElim}(i) \cup {i} = Adj_{GElim}(j) \cup {j},
  //
  // but we instead use the stricter, and easier to measure, test that
  //
  //   Adj_{GQuotient}(i) \cup {i} = Adj_{GQuotient}(j) \cup {j}.
  //
  // The second test is easier to compute because
  //
  //   Adj_{GQuotient}(i) \cup {i} = A_i \cup E_i \cup {i},
  //
  // and
  //
  //   A_i \cup E_i \cup {i} = A_j \cup E_j \cup {j}
  //
  // if and only if
  //
  //   A_i \cup {i} = A_j \cup {j} and E_i = E_j.
  //
  bool VariablesAreQuotientIndistinguishable(Int i, Int j) const;

  // An implementation of Algorithm 2 from [ADD-96].
  // On exit, it holds |L_e \ L_p| for all elements e in the element list
  // of a supernode in the structure, L_p.
  std::unordered_map<Int, Int> ExternalStructureSizes(
    const std::vector<Int>& supernodal_structure) const;
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
      result = (result + index) % (num_original_vertices - 1);
    }
  }
  for (const Int& index : element_lists[i]) {
    result = (result + index) % (num_original_vertices - 1);
  }
  return result + 1;
};

inline bool QuotientGraph::VariablesAreQuotientIndistinguishable(
    Int i, Int j) const {
  const std::size_t adj_size = adjacency_lists[i].size();
  const std::size_t elem_size = element_lists[i].size();

  // Early exit if the set cardinalities disagree.
  if (adj_size != adjacency_lists[j].size() ||
      elem_size != element_lists[j].size()) {
    return false;
  }

  // Check if E_i = E_j.
  for (std::size_t index = 0; index < elem_size; ++index) {
    if (element_lists[i][index] != element_lists[j][index]) {
      return false;
    }
  }

  // Check if A_i \cup {i} = A_j \cup {j}.
  const std::vector<Int> i_vec{i}, j_vec{j};
  std::vector<std::vector<Int>::const_iterator> i_iters{
    adjacency_lists[i].cbegin(), i_vec.cbegin(),
  };
  std::vector<std::vector<Int>::const_iterator> j_iters{
    adjacency_lists[j].cbegin(), j_vec.cbegin(),
  };
  const std::vector<std::vector<Int>::const_iterator> i_ends{
    adjacency_lists[i].cend(), i_vec.cend(),
  };
  const std::vector<std::vector<Int>::const_iterator> j_ends{
    adjacency_lists[j].cend(), j_vec.cend(),
  };

  for (std::size_t index = 0; index < adj_size + 1; ++index) {
    bool found_match = false;
    for (std::size_t i_set = 0; i_set < i_iters.size(); ++i_set) {
      auto& i_iter = i_iters[i_set];
      if (i_iter == i_ends[i_set]) {
        continue;
      }

      for (std::size_t j_set = 0; j_set < j_iters.size(); ++j_set) {
        auto& j_iter = j_iters[j_set];
        if (j_iter == j_ends[j_set]) {
          continue;
        }

        if (*i_iter == *j_iter) {
          ++i_iter;
          ++j_iter;
          found_match = true;
          break;
        }
      }
      if (found_match) {
        break;
      }
    }
    if (!found_match) {
      return false;
    }
  }
#ifdef QUOTIENT_DEBUG
  for (Int i_set = 0; i_set < 2; ++i_set) {
    if (i_iters[i_set] != i_ends[i_set]) {
      std::cerr << "i_iters[" << i_set << "] != end." << std::endl;
    }
  }
  for (Int j_set = 0; j_set < 2; ++j_set) {
    if (j_iters[j_set] != j_ends[j_set]) {
      std::cerr << "j_iters[" << j_set << "] != end." << std::endl;
    }
  }
#endif

  return true;
}

inline std::unordered_map<Int, Int> QuotientGraph::ExternalStructureSizes(
    const std::vector<Int>& supernodal_structure) const {
  std::unordered_map<Int, Int> external_structure_sizes;
  for (const Int& i : supernodal_structure) {
    for (const Int& element : element_lists[i]) {
      if (!external_structure_sizes.count(element)) {
        external_structure_sizes[element] = structures[element].size();
      }
      external_structure_sizes[element] -= supernodes[i].size();
    }
  }
  return external_structure_sizes;
}

} // namespace quotient

#endif // ifndef QUOTIENT_QUOTIENT_GRAPH_H_
