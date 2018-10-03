/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_MINIMUM_DEGREE_H_
#define QUOTIENT_MINIMUM_DEGREE_H_

#include <iostream>
#include <unordered_map>
#include <unordered_set>
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
// There appears to be a plethora of publicly available preprints. We 
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

// A sequence of external degree approximations of decreasing accuracy.
// Please see Theorem 4.1 of [ADD-96] for accuracy guarantees.
enum ExternalDegreeType {
  // In the notation of [ADD-96],
  //   d_i = |A_i \ supernode(i)| + |(\cup_{e in E_i} L_e) \ supernode(i)|.
  kExactExternalDegree,

  // In the notation of [ADD-96],
  //   \bar{d_i}^k = min(
  //     n - k,
  //     \bar{d_i}^{k - 1} + |L_p \ supernode(i)|, 
  //     |A_i \ supernode(i)| + |L_p \ supernode(i)| +
  //         \sum_{e in E_i \ {p}} |L_e | L_p|).
  kAmestoyExternalDegree,

  // In the notation of [ADD-96],
  //   \tilde{d_i} = d_i if |E_i| = 2, \hat{d_i}, otherwise.
  kAshcraftExternalDegree,

  // In the notation of [ADD-96],
  //   \hat{d_i} = |A_i \ supernode(i)| + \sum_{e in E_i} |L_e \ supernode(i)|.
  kGilbertExternalDegree,
};

// Returns the number of entries in {vec0} \cap {vec1}, where both vectors
// are assumed sorted and unique.
template<typename T>
Int SizeOfIntersection(
    const std::vector<T>& vec0, const std::vector<T>& vec1) {
  Int num_intersections = 0;
  auto vec0_iter = vec0.cbegin(); 
  auto vec1_iter = vec1.cbegin();
  while (vec0_iter != vec0.cend() && vec1_iter != vec1.cend()) {
    if (*vec0_iter < *vec1_iter) {
      vec0_iter = std::lower_bound(vec0_iter, vec0.cend(), *vec1_iter);
    } else if (*vec0_iter > *vec1_iter) {
      vec1_iter = std::lower_bound(vec1_iter, vec1.cend(), *vec0_iter);
    } else {
      ++num_intersections;
      ++vec0_iter;
      ++vec1_iter;
    }
  }
  return num_intersections;
}

// Returns the number of entries in {vec} \ {blacklist}, where both vectors are
// assumed to be sorted and unique.
template<typename T>
Int SizeOfDifference(
    const std::vector<T>& vec, const std::vector<T>& blacklist) {
  return vec.size() - SizeOfIntersection(vec, blacklist);
}

// Returns the number of entries in {vec0} \cup {vec1}, where both vectors are
// assumed to be sorted and unique.
template<typename T>
Int SizeOfUnion(const std::vector<T>& vec0, const std::vector<T>& vec1) {
  return vec0.size() + vec1.size() - SizeOfIntersection(vec0, vec1);
}

// Returns the number of entries in ({vec0} \cup {vec1}) \ {blacklist}, where
// all vectors are assumed sorted and unique.
template<typename T>
Int SizeOfBlacklistedUnion(
    const std::vector<T>& vec0,
    const std::vector<T>& vec1,
    const std::vector<T>& blacklist) {
  auto vec0_iter = vec0.cbegin(); 
  auto blacklist_iter0 = blacklist.cbegin();

  auto vec1_iter = vec1.cbegin();
  auto blacklist_iter1 = blacklist.cbegin();

  Int num_vec0_blacklisted = 0;
  Int num_vec1_blacklisted = 0;
  Int num_blacklisted_intersections = 0;
  while (vec0_iter != vec0.cend() && vec1_iter != vec1.cend()) {
    // Skip this entry of vec0 if it is blacklisted.
    blacklist_iter0 =
      std::lower_bound(blacklist_iter0, blacklist.cend(), *vec0_iter);
    if (blacklist_iter0 != blacklist.cend() && *vec0_iter == *blacklist_iter0) {
      ++vec0_iter; 
      ++blacklist_iter0;
      ++num_vec0_blacklisted;
      continue;
    }

    // Skip this entry of vec1 if it is blacklisted.
    blacklist_iter1 =
      std::lower_bound(blacklist_iter1, blacklist.cend(), *vec1_iter);
    if (blacklist_iter1 != blacklist.cend() && *vec1_iter == *blacklist_iter1) {
      ++vec1_iter;
      ++blacklist_iter1;
      ++num_vec1_blacklisted;
      continue;
    }

    if (*vec0_iter < *vec1_iter) {
      vec0_iter = std::lower_bound(vec0_iter, vec0.cend(), *vec1_iter);
    } else if (*vec0_iter > *vec1_iter) {
      vec1_iter = std::lower_bound(vec1_iter, vec1.cend(), *vec0_iter);
    } else {
      ++num_blacklisted_intersections;
      ++vec0_iter;
      ++vec1_iter;
    }
  }

  return (vec0.size() - num_vec0_blacklisted) +
      (vec1.size() - num_vec1_blacklisted) - num_blacklisted_intersections;
}

// Fills 'filtered_vec' with the sorted serialization of {vec} \ {blacklist},
// where 'vec' and 'blacklist' are both sorted, unique vectors.
template<typename T>
void FilterSet(
    const std::vector<T>& vec,
    const std::vector<T>& blacklist,
    std::vector<T>* filtered_vec) {
  const Int filtered_size = SizeOfDifference(vec, blacklist);
  filtered_vec->resize(0);
  filtered_vec->resize(filtered_size);
#ifdef QUOTIENT_DEBUG
  auto iter = std::set_difference(
      vec.begin(), vec.end(),
      blacklist.begin(), blacklist.end(),
      filtered_vec->begin());
  if (filtered_size != std::distance(filtered_vec->begin(), iter)) {
    std::cerr << "Filtered sizes did not match in FilterSet." << std::endl;
  }
#else
  std::set_difference(
      vec.begin(), vec.end(),
      blacklist.begin(), blacklist.end(),
      filtered_vec->begin());
#endif
}

// Fills 'sorted_union' with the sorted serialization of {vec0} \cup {vec1},
// where 'vec0' and 'vec1' are sorted, unique vectors.
template<typename T>
void MergeSets(
    const std::vector<T>& vec0,
    const std::vector<T>& vec1, 
    std::vector<T>* sorted_union) {
  const Int union_size = SizeOfUnion(vec0, vec1);
  sorted_union->resize(0);
  sorted_union->resize(union_size);
#ifdef QUOTIENT_DEBUG
  auto iter = std::set_union(
      vec0.begin(), vec0.end(),
      vec1.begin(), vec1.end(),
      sorted_union->begin());
  if (union_size != std::distance(sorted_union->begin(), iter)) {
    std::cerr << "Union sizes did not match in MergeSets." << std::endl;
  }
#else
  std::set_union(
      vec0.begin(), vec0.end(),
      vec1.begin(), vec1.end(),
      sorted_union->begin());
#endif
}

// Inserts 'value' into the sorted, unique vector, 'vec', so that the result
// is sorted.
template<typename T>
void InsertEntryIntoSet(const T& value, std::vector<T>* vec) {
  auto iter = std::lower_bound(vec->begin(), vec->end(), value);
  vec->insert(iter, value);
}

// A definition of Ashcraft's hash function (as described in [ADD-96]).
// Note that only principal members of A_i are incorporated in the hash.
inline std::size_t AshcraftVariableHash(const QuotientGraph& graph, Int i) {
  std::size_t result = 0;
  for (const Int& index : graph.adjacency_lists[i]) {
    if (!graph.supernodes[index].empty()) {
      result = (result + index) % (graph.num_original_vertices - 1);
    }
  }
  for (const Int& index : graph.element_lists[i]) {
    result = (result + index) % (graph.num_original_vertices - 1);
  }
  return result + 1;
};

// Returns true if supernodes 'i' and 'j' are considered indistinguishable with
// respect to their quotient graph representation.
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
inline bool SupernodesAreQuotientIndistinguishable(
    const QuotientGraph& graph, Int i, Int j) {
  const std::size_t adj_size = graph.adjacency_lists[i].size();
  const std::size_t elem_size = graph.element_lists[i].size();

  // Early exit if the set cardinalities disagree.
  if (adj_size != graph.adjacency_lists[j].size() ||
      elem_size != graph.element_lists[j].size()) {
    return false;
  }

  // Check if E_i = E_j.
  for (std::size_t index = 0; index < elem_size; ++index) {
    if (graph.element_lists[i][index] != graph.element_lists[j][index]) {
      return false;
    }
  }

  // Check if A_i \cup {i} = A_j \cup {j}.
  const std::vector<Int> i_vec{i}, j_vec{j};
  std::vector<std::vector<Int>::const_iterator> i_iters{
    graph.adjacency_lists[i].cbegin(), i_vec.cbegin(),
  };
  std::vector<std::vector<Int>::const_iterator> j_iters{
    graph.adjacency_lists[j].cbegin(), j_vec.cbegin(),
  };
  const std::vector<std::vector<Int>::const_iterator> i_ends{
    graph.adjacency_lists[i].cend(), i_vec.cend(),
  };
  const std::vector<std::vector<Int>::const_iterator> j_ends{
    graph.adjacency_lists[j].cend(), j_vec.cend(),
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

// An implementation of Algorithm 2 from [ADD-96].
// On exit, it holds |L_e \ L_p| for all elements e in the element list
// of a supernode in the pivot structure.
inline std::unordered_map<Int, Int> ExternalStructureSizes(
    const QuotientGraph& graph,
    const std::vector<Int>& supernodal_pivot_structure) {
  std::unordered_map<Int, Int> external_structure_sizes;
  for (const Int& i : supernodal_pivot_structure) {
    for (const Int& element : graph.element_lists[i]) {
      if (!external_structure_sizes.count(element)) {
        external_structure_sizes[element] = graph.structures[element].size();
      }
      external_structure_sizes[element] -= graph.structures[i].size();
    }
  }
  return external_structure_sizes;
}

// Computes the exact external degree of supernode i using Eq. (2) of [ADD-96].
//   d_i = |A_i \ supernode(i)| + |(\cup_{e in E_i) L_e) \ supernode(i)|.
inline Int ExactExternalDegree(const QuotientGraph& graph, Int i) {
  // Add the cardinality of A_i \ supernode(i).
  Int external_degree = SizeOfDifference(
      graph.adjacency_lists[i], graph.supernodes[i]);

  // Add the cardinality of (\cup_{e in E_i} L_e) \ supernode(i).
  if (graph.element_lists[i].size() == 2) {
    const Int element0 = graph.element_lists[i][0];
    const Int element1 = graph.element_lists[i][1];
    external_degree += SizeOfBlacklistedUnion(
      graph.structures[element0],
      graph.structures[element1],
      graph.supernodes[i]);
  } else {
    std::vector<Int> filtered_struct;
    std::vector<Int> element_struct_union;
    std::vector<Int> temp_int_vec;
    for (const Int& element : graph.element_lists[i]) {
      FilterSet(graph.structures[element], graph.supernodes[i],
                &filtered_struct);
      auto temp_int_vec = element_struct_union;
      MergeSets(temp_int_vec, filtered_struct, &element_struct_union);
    }
    external_degree += element_struct_union.size();
  }

  return external_degree;
}

// Computes an approximation of the external degree of supernode i using Eq. (4)
// of [ADD-96].
inline Int AmestoyExternalDegree(
    const QuotientGraph& graph,
    Int i,
    Int pivot,
    Int old_external_degree,
    const std::unordered_map<Int, Int>& external_structure_sizes) {
  const Int num_vertices_left =
      graph.num_original_vertices - graph.num_eliminated_vertices;

  // Note that this usage of 'external' refers to |L_p \ supernode(i)| and not
  // |L_e \ L_p|, as is the case for 'external_structure_sizes'.
  const Int external_pivot_structure_size = SizeOfDifference(
      graph.structures[pivot], graph.supernodes[i]);

  const Int bound0 = num_vertices_left;
  const Int bound1 = old_external_degree + external_pivot_structure_size;

  Int bound2 =
    SizeOfDifference(graph.adjacency_lists[i], graph.supernodes[i]) +
    SizeOfDifference(graph.structures[pivot], graph.supernodes[i]);
  for (const Int& element : graph.element_lists[i]) {
    if (element == pivot) {
      continue;
    }
    if (external_structure_sizes.count(element)) {
      bound2 += graph.structures[element].size();
    } else {
      auto iter = external_structure_sizes.find(element);
#ifdef QUOTIENT_DEBUG
      if (iter == external_structure_sizes.end()) {
        std::cerr << "Did not find element." << std::endl;
      }
#endif
      bound2 += iter->second;
    }
  }

  return std::min(bound0, std::min(bound1, bound2));
}

// Returns the external degree approximation of Gilbert, Moler, and Schreiber,
//   \hat{d_i} = |A_i \ supernode(i)|  + \sum_{e in E_i} |L_e \ supernode(i)|.
inline Int GilbertExternalDegree(const QuotientGraph& graph, Int i) {
  Int external_degree = SizeOfDifference(
      graph.adjacency_lists[i], graph.supernodes[i]);
  for (const Int& element : graph.element_lists[i]) {
    external_degree += SizeOfDifference(
        graph.structures[element], graph.supernodes[i]);
  }
  return external_degree;
}

// Returns the external degree approximation of Ashcraft, Eisenstat, and Lucas:
//   \tilde{d_i} = d_i if |E_i| = 2, \hat{d_i} otherwise.
inline Int AshcraftExternalDegree(const QuotientGraph& graph, Int i) {
  if (graph.element_lists[i].size() == 2) {
    return ExactExternalDegree(graph, i);
  }
  return GilbertExternalDegree(graph, i);
}

// Returns (an approximation of) the external degree of a given supervariable.
inline Int ExternalDegree(
    const QuotientGraph& graph,
    Int i,
    Int pivot,
    Int old_external_degree,
    const std::unordered_map<Int, Int>& external_structure_sizes,
    ExternalDegreeType degree_type) {
  Int external_degree = -1;
  switch(degree_type) {
    case kExactExternalDegree: { 
      external_degree = ExactExternalDegree(graph, i);
      break;
    }
    case kAmestoyExternalDegree: {
      external_degree = AmestoyExternalDegree(
          graph, i, pivot, old_external_degree, external_structure_sizes);
      break;
    }
    case kAshcraftExternalDegree: { 
      external_degree = AshcraftExternalDegree(graph, i);
      break;
    }
    case kGilbertExternalDegree: {
      external_degree = GilbertExternalDegree(graph, i);
      break;
    }
  }
  return external_degree;
}

// The result of running the MinimumDegree reordering algorithm. It contains
// the ordered list of eliminated principal vertices, the list of supernodes,
// and the supernodal nonzero structure of each principal column.
struct MinimumDegreeAnalysis {
  // The recommended elimination order of the supernodes (with each supernode
  // represented by its principal member).
  std::vector<Int> elimination_order;

  // The list of supernodes. Entry 'i' is nonempty if and only if 'i' is the
  // principal member of the supernode, in which case 'supernodes[i]' contains
  // the sorted list of members of the supernode.
  std::vector<std::vector<Int>> supernodes;

  // The supernodal nonzero structure of the principal columns of the
  // factorization.
  std::vector<std::vector<Int>> supernodal_structures;

  // We will push to elimination order as the reordering algorithm progresses,
  // so we will allocate an upper bound for the amount of required space.
  // The 'supernodes' and 'structures' variables will be copied over from
  // the quotient graph just before the analysis completes.
  MinimumDegreeAnalysis(Int num_vertices) {
    elimination_order.reserve(num_vertices);
  }
};

// Returns a supernodal reordering and the corresponding supernodal nonzero
// structure of the implied factorization using the (Approximate) Minimum
// Degree reordering algorithm. Please see [ADD-96] for details.
//
// The input graph must be explicitly symmetric.
//
MinimumDegreeAnalysis MinimumDegree(
  const CoordinateGraph& graph, ExternalDegreeType degree_type) {
#ifdef QUOTIENT_DEBUG
  if (graph.NumSources() != graph.NumVertices()) {
    std::cerr << "ERROR: MinimumDegree requires a symmetric input graph."
              << std::endl;
    return MinimumDegreeAnalysis(0);
  }
#endif
  QuotientGraph quotient_graph(graph);
  const Int num_orig_vertices = quotient_graph.num_original_vertices;

  // Initialize the cached binary tree of external degrees.
  RandomAccessHeap<Int> external_degree_heap;
  {
    std::vector<Int> external_degrees_vec(num_orig_vertices);
    for (Int source = 0; source < num_orig_vertices; ++source) {
      external_degrees_vec[source] =
          quotient_graph.adjacency_lists[source].size();
    }
    external_degree_heap.Reset(external_degrees_vec);
  }

  // Buffers used to store a temporary integer vector.
  // They are declared here to avoid unnecessary garbage collection.
  std::vector<Int> temp_int_vec0, temp_int_vec1;

  // Eliminate the variables.
  MinimumDegreeAnalysis analysis(num_orig_vertices);
  while (quotient_graph.num_eliminated_vertices < num_orig_vertices) {
    // Retrieve a variable with minimal (approximate) external degree.
    const std::pair<Int, Int> pivot_pair = external_degree_heap.MinimalEntry(); 
    const Int pivot = pivot_pair.first;

    // Compute the structure of the pivot:
    //   L_p := (A_p \cup (\cup_{e in E_p} L_e)) \ p.
    //
    // TODO(Jack Poulson): Experiment with several different union-finding
    // algorithms. It isn't clear that repeated calls to std::set_union has
    // optimal time complexity.
    FilterSet(
      quotient_graph.adjacency_lists[pivot],
      quotient_graph.supernodes[pivot],
      &quotient_graph.structures[pivot]);
    for (const Int& element : quotient_graph.element_lists[pivot]) {
      FilterSet(
        quotient_graph.structures[element],
        quotient_graph.supernodes[pivot],
        &temp_int_vec0);

      temp_int_vec1 = quotient_graph.structures[pivot];
      MergeSets(
        temp_int_vec1,
        temp_int_vec0,
        &quotient_graph.structures[pivot]);
    }
    const std::vector<Int> supernodal_pivot_structure =
        quotient_graph.FormSupernodalStructure(pivot);

    // Update the adjacency lists.
    for (const Int& i : supernodal_pivot_structure) {
      // Remove redundant adjacency entries:
      //   A_i := (A_i \ L_p) \ supernode(p).
      FilterSet(
          quotient_graph.adjacency_lists[i],
          quotient_graph.structures[pivot],
          &temp_int_vec0);
      FilterSet(
          temp_int_vec0,
          quotient_graph.supernodes[pivot],
          &quotient_graph.adjacency_lists[i]);
    }

    // Update the element lists.
    for (const Int& i : supernodal_pivot_structure) {
      // Element absorption:
      //   E_i := (E_i \ E_p) \cup {p}.
      temp_int_vec0 = quotient_graph.element_lists[i];
      FilterSet(
          temp_int_vec0,
          quotient_graph.element_lists[pivot],
          &quotient_graph.element_lists[i]);
      InsertEntryIntoSet(pivot, &quotient_graph.element_lists[i]);
    }

    // Compute the external structure cardinalities of the elements.
    // (but only if the Amestoy external degree approximation is requested).
    std::unordered_map<Int, Int> external_structure_sizes;
    if (degree_type == kAmestoyExternalDegree) {
      external_structure_sizes = ExternalStructureSizes(
          quotient_graph, supernodal_pivot_structure);
    }

    // Compute the external degree approximations of the supernodes
    // adjacent to the current pivot.
    //
    // TODO(Jack Poulson): Add support for a batch interface to
    // 'external_degree_heap.SetValue' so that comparison propagations can be
    // potentially shared.
    for (const Int& i : supernodal_pivot_structure) {
      // Compute the external degree (or approximation) of supervariable i:
      //   d_i := |A_i \ supernode(i)| + |(\cup_{e in E_i} L_e) \ supernode(i)|.
      const Int old_external_degree = external_degree_heap.Value(i); 
      const Int external_degree = ExternalDegree(
        quotient_graph, i, pivot, old_external_degree, external_structure_sizes,
        degree_type);
      external_degree_heap.SetValue(i, external_degree);
    }

    // Fill a set of buckets for the hashes of the supernodes adjacent to
    // the current pivot.
    const Int supernodal_pivot_struct_size = supernodal_pivot_structure.size();
    std::vector<Int> bucket_list(supernodal_pivot_struct_size);
    std::unordered_map<Int, std::vector<Int>> variable_hash_map;
    for (Int i_index = 0; i_index < supernodal_pivot_struct_size; ++i_index) {
      const Int i = supernodal_pivot_structure[i_index];

      // Append this principal variable to its hash bucket.
      // TODO(Jack Poulson): Add configurable support for other hashes.
      // For example, one could mod by (std::numeric_limits<Int>::max() - 1).
      const std::size_t bucket = AshcraftVariableHash(quotient_graph, i);
      bucket_list[i_index] = bucket;
      variable_hash_map[bucket].push_back(i_index);
    }

    // Supervariable detection. While the supernodal merges will potentially
    // shrink the supernodal adjacency lists (and thus change the associated
    // Ashcraft hash of variables), if two variables are indistinguishable,
    // their cached bucket might be wrong, but they would be wrong together.
    //
    // The test for indistinguishability does not depend upon the variable
    // supernodal structure and is thus invariant to supervariable merges.
    std::vector<bool> merged_supernode(supernodal_pivot_struct_size, false);
    for (Int i_index = 0; i_index < supernodal_pivot_struct_size; ++i_index) {
      if (merged_supernode[i_index]) {
        continue;
      }
      const Int i = supernodal_pivot_structure[i_index];
      const Int bucket = bucket_list[i_index];
      for (const Int& j_index : variable_hash_map[bucket]) {
        if (j_index <= i_index || merged_supernode[j_index]) {
          continue;
        }
        const Int j = supernodal_pivot_structure[j_index];
        const bool indistinguishable =
            SupernodesAreQuotientIndistinguishable(quotient_graph, i, j);
        if (indistinguishable) {
          // Absorb supernode(j) into supernode(i). 
          temp_int_vec0 = quotient_graph.supernodes[i];
          MergeSets(
              temp_int_vec0,
              quotient_graph.supernodes[j],
              &quotient_graph.supernodes[i]);
          external_degree_heap.UpdateValue(
              i, -quotient_graph.supernodes[j].size());
          external_degree_heap.DisableIndex(j);
#ifdef QUOTIENT_DEBUG
          quotient_graph.variables.erase(j);
#endif
          SwapClearVector(&quotient_graph.supernodes[j]);
          SwapClearVector(&quotient_graph.adjacency_lists[j]);
          SwapClearVector(&quotient_graph.element_lists[j]);
          merged_supernode[j_index] = true;
        }
      }
    }

    // Convert pivot from a variable to an element.
#ifdef QUOTIENT_DEBUG
    // Update the element list:
    //   \bar{V} := (\bar{V} \cup {p}) \ E_p
    quotient_graph.elements.insert(pivot);
    for (const Int& element : quotient_graph.element_lists[pivot]) {
      quotient_graph.elements.erase(element);
    }
    // Update the variable list:
    //   V := V \ {p}.
    quotient_graph.variables.erase(pivot);
#endif
    external_degree_heap.DisableIndex(pivot);
    SwapClearVector(&quotient_graph.adjacency_lists[pivot]);
    SwapClearVector(&quotient_graph.element_lists[pivot]);
    quotient_graph.num_eliminated_vertices +=
        quotient_graph.supernodes[pivot].size();

    analysis.elimination_order.push_back(pivot);
  }

  analysis.supernodes = quotient_graph.supernodes;
  analysis.supernodal_structures.resize(num_orig_vertices);
  for (const Int& i : analysis.elimination_order) {
    analysis.supernodal_structures[i] =
        quotient_graph.FormSupernodalStructure(i);
  }
  return analysis;
}

} // namespace quotient

#endif // ifndef QUOTIENT_MINIMUM_DEGREE_H_
