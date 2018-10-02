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
#include "quotient/shrinkable_cached_binary_tree.hpp"

namespace quotient {

// TODO(Jack Poulson): Carefully describe this data structure.
struct QuotientGraph {
  Int num_orignal_vertices;

  Int num_eliminated_vertices;

  std::vector<std::vector<Int>> supervariables;

  std::vector<std::vector<Int>> structures;

  std::vector<std::vector<Int>> adjacency_lists;

  std::vector<std::vector<Int>> element_lists;

#ifdef QUOTIENT_DEBUG
  std::set<Int> variables;

  std::set<Int> elements;
#endif

  QuotientGraph(const CoordinateGraph& graph)
  : num_original_vertices(graph.NumSources()), num_eliminated_vertices(0) {
    // Initialize the supervariables as simple.
    supervariables.resize(num_original_vertices);
    for (Int source = 0; source < num_original_vertices; ++source) {
      supervariables[source].push_back(source);
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
};

// A sequence of external degree approximations of decreasing accuracy.
enum ExternalDegreeType {
  // In the notation of Amestoy et al. 1996,
  //   d_i = |A_i \ supervar(i)| + |(\cup_{e in E_i} L_e) \ supervar(i)|.
  kExactExternalDegree,

  // In the notation of Amestoy et al. 1996,
  //   \bar{d_i}^k = min(
  //     n - k,
  //     \bar{d_i}^{k - 1} + |L_p \ supervar(i)|, 
  //     |A_i \ supervar(i)| + |L_p \ supervar(i)| +
  //         \sum_{e in E_i \ {p}} |L_e | L_p|).
  kAmestoyExternalDegree,

  // In the notation of Amestoy et al. 1996,
  //   \tilde{d_i} = d_i if |E_i| = 2, \hat{d_i}, otherwise.
  kAshcraftExternalDegree,

  // In the notation of Amestoy et al. 1996,
  //   \hat{d_i} = |A_i \ supervar(i)| + \sum_{e in E_i} |L_e \ supervar(i)|.
  kGilbertExternalDegree,
};

// Returns the number of entries in {vec} \cap {blacklist}, where both vectors
// are assumed sorted and unique.
template<typename T>
Int CardinalityOfIntersection(
    const std::vector<T>& vec, const std::vector<T>& blacklist) {
  Int num_intersections = 0;
  auto vec_iter = vec.cbegin(); 
  auto blacklist_iter = blacklist.cbegin();
  while (vec_iter != vec.cend() && blacklist_iter != blacklist.cend()) {
    if (*vec_iter < *blacklist_iter) {
      vec_iter = std::lower_bound(vec_iter, vec.cend(), *blacklist_iter);
    } else if (*vec_iter > *blacklist_iter) {
      blacklist_iter = std::lower_bound(
          blacklist_iter, blacklist.cend(), *vec_iter);
    } else {
      ++num_removed;
      ++vec_iter;
      ++blacklist_iter;
    }
  }
  return num_intersections;
}

// Returns the sorted serialization of {vec} \ {blacklist}, where 'vec' and
// 'blacklist' are both sorted, unique vectors.
template<typename T>
void FilterSet(
    const std::vector<T>& vec,
    const std::vector<T>& blacklist,
    std::vector<T>* filtered_vec) {
  const Int num_intersections = CardinalityOfIntersection(vec, blacklist);
  const Int filtered_size = vec.size() - num_intersections;
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

template<typename T>
void MergeSets(
    const std::vector<T>& vec0,
    const std::vector<T>& vec1, 
    std::vector<T>* sorted_union) {
  const Int num_intersections = CardinalityOfIntersection(vec0, vec1);
  const Int union_size = vec0.size() + vec1.size() - num_intersections;
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

inline void InsertEntryIntoSet(Int value, std::vector<Int>* vec) {
  auto iter = std::lower_bound(vec->begin(), vec->end(), value);
  vec->insert(iter, value);
}

// A definition of Ashcraft's hash function (as described in Amestoy et al.,
// "An Approximate Minimum Degree Ordering Algorithm").
inline std::size_t AshcraftVariableHash(
    Int num_orig_vertices,
    const std::vector<Int>& adjacency_list,
    const std::vector<Int>& element_list) {
  std::size_t result = 0;
  for (const Int& index : adjacency_list) {
    result = (result + index) % (num_orig_vertices - 1);
  }
  for (const Int& index : element_list) {
    result = (result + index) % (num_orig_vertices - 1);
  }
  return result + 1;
};

// Returns true if supernodes 'i' and 'i' are considered indistinguishable with
// respect to their quotient graph representation.
//
// The elimination graph definition (e.g., as given by Amestoy et al.)
// involves testing if
//   Adj_{GElim}(i) \cup {i} = Adj_{GElim}(j) \cup {j},
// but we instead use the stricter, and easier to measure, test that
//   Adj_{GQuotient}(i) \cup {i} = Adj_{GQuotient}(j) \cup {j}.
//
// The second test is easier to compute because
///  Adj_{GQuotient}(i) \cup {i} = A_i \cup E_i \cup {i},
// where A_i, E_i, and {i} are nonoverlapping, sorted lists.
//
// We can therefore test for indistinguishability by advancing three pointers
// for each supervariable.
inline bool SupernodesAreQuotientIndistinguishable(
    const QuotientGraph& graph, Int i, Int j) {
  const Int adj_i_size = graph.adjacency_lists[i].size();
  const Int adj_j_size = graph.adjacency_lists[j].size();
  const Int elem_i_size = graph.element_lists[i].size();
  const Int elem_j_size = graph.element_lists[j].size();
  const Int set_size = adj_i_size + elem_i_size + 1;

  // Early exit if the set cardinalities disagree.
  if (set_size != adj_j_size + elem_j_size + 1) {
    return false;
  }

  // Set up iterator bounds for the adjacency list, element list, and singleton.
  const std::vector<Int> i_vec{i};
  const std::vector<Int> j_vec{j}; 
  std::vector<std::vector<Int>::const_iterator> i_iters{
    graph.adjacency_lists[i].cbegin(),
    graph.element_lists[i].cbegin(),
    i_vec.cbegin(),
  };
  std::vector<std::vector<Int>::const_iterator> j_iters{
    graph.adjacency_lists[j].cbegin(),
    graph.element_lists[j].cbegin(),
    j_vec.cbegin(),
  };
  const std::vector<std::vector<Int>::const_iterator> i_ends{
    graph.adjacency_lists[i].cend(),
    graph.element_lists[i].cend(),
    i_vec.cend(),
  };
  const std::vector<std::vector<Int>::const_iterator> j_ends{
    graph.adjacency_lists[j].cend(),
    graph.element_lists[j].cend(),
    j_vec.cend(),
  };
 
  for (Int index = 0; index < set_size; ++index) {
    bool found_match = false;
    for (Int i_set = 0; i_set < 3; ++i_set) {
      if (i_iters[i_set] != i_ends[i_set]) {
        for (Int j_set = 0; j_set < 3; ++j_set) {
          if (j_iters[j_set] != j_ends[j_set] &&
              *i_iters[i_set] == *j_iters[j_set]) {
            ++i_iters[i_set];
            ++j_iters[j_set];
            found_match = true;
            break;
          }
        }
        if (found_match) {
          break;
        }
      }
    }
    if (!found_match) {
      return false;
    }
  }
#ifdef QUOTIENT_DEBUG
  for (Int i_set = 0; i_set < 3; ++i_set) {
    if (i_iters[i_set] != i_ends[i_set]) {
      std::cerr << "i_iters[" << i_set << "] != end." << std::endl;
    }
  }
  for (Int j_set = 0; j_set < 3; ++j_set) {
    if (j_iters[j_set] != j_ends[j_set]) {
      std::cerr << "j_iters[" << j_set << "] != end." << std::endl;
    }
  }
#endif

  return true;
}

// An implementation of Algorithm 2 from Amestoy et al., 1996.
// It implicitly computes |L_e \ L_p| for all elements e within an
// unordered map from the element index to the cardinality. If no key exists
// for a particular element index, the cardinality is implied to be |L_e|.
inline std::unordered_map<Int, Int> ExternalStructureCardinalities(
    const QuotientGraph& graph, Int pivot) {
  std::unordered_map<Int, Int> external_structure_cardinalities;
  for (const Int& i : graph.structures[pivot]) {
    for (const Int& element : graph.element_lists[i]) {
      if (!external_structure_cardinalities.count(element)) {
        external_structure_cardinalities[element] =
            graph.structures[element].size();
      }
      external_structure_cardinalities[element] -= graph.structures[i].size();
    }
  }
  return external_structure_cardinalities;
}

inline Int ExactExternalDegree(const QuotientGraph& graph, Int i) {
  // TODO(Jack Poulson)
}

inline Int AmestoyExternalDegree(
    const QuotientGraph& graph,
    Int i,
    Int pivot,
    Int old_external_degree,
    const std::unordered_map<Int, Int>& external_structure_cardinalities) {
  const Int num_vertices_left =
      graph.num_original_vertices - graph.num_eliminated_vertices;
  // TODO(Jack Poulson)
}

// Returns the external degree approximation of Gilbert, Moler, and Schreiber,
//   \hat{d_i} = |A_i \ supervar(i)|  + \sum_{e in E_i} |L_e \ supervar(i)|.
inline Int GilbertExternalDegree(const QuotientGraph& graph, Int i) {
  Int external_degree = graph.adjacency_lists[i].size() -
      CardinalityOfIntersection(
          graph.adjacency_lists[i],
          graph.supervariables[i]);
  for (const Int& element : element_list) {
    external_degree += graph.structures[element] -
        CardinalityOfIntersection(
            graph.structures[element],
            graph.supervariables[i]);
  }
  return external_degree;
}

// Returns the external degree approximation of Ashcraft, Eisenstat, and Lucas:
//   \tilde{d_i} = d_i if |E_i| = 2, \hat{d_i} otherwise.
inline Int AshcraftExternalDegree(const QuotientGraph& graph, Int i) {
  if (graph.element_lists[i].size() == 2) {
    // TODO(Jack Poulson)
  }
  return GilbertExternalDegree(graph, i);
}

// Returns (an approximation of) the external degree of a given supervariable.
inline Int ExternalDegree(
    const QuotientGraph& graph,
    Int i,
    Int pivot,
    Int old_external_degree,
    const std::unordered_map<Int, Int>& external_structure_cardinalities,
    ExternalDegreeType degree_type) {
  switch(degree_type) {
    case kExactExternalDegree: { 
      return ExactExternalDegree(graph, i);
    }
    case kAmestoyExternalDegree: {
      return AmestoyExternalDegree(
          graph, i, pivot, old_external_degree,
          external_structure_cardinalities);
    }
    case kAshcraftExternalDegree: { 
      return AshcraftExternalDegree(graph, i);
    }
    case kGilbertExternalDegree: {
      return GilbertExternalDegree(graph, i);
    }
  }
}

struct MinimumDegreeAnalysis {
  std::vector<Int> elimination_order;

  std::vector<std::vector<Int>> supervariables;

  std::vector<std::vector<Int>> structures;

  // We will push to elimination order as the reordering algorithm progresses,
  // so we will allocate an upper bound for the amount of required space.
  // The 'supervariables' and 'structures' variables will be copied over from
  // the quotient graph just before the analysis completes.
  Analysis(Int num_vertices) {
    elimination_order.reserve(num_vertices);
  }
};

// Returns a list of supernodes of the given (symmetric) graph in an order
// meant to approximately minimize fill-in during Cholesky factorization.
// The input graph must be explicitly symmetric.
//
// Please see:
//   Patrick R. Amestoy, Timothy A. Davis, and Iain S. Duff,
//   "An Approximate Minimum Degree Ordering Algorithm",
//   SIAM J. Matrix Analysis & Applic., Vol. 17, No. 4, pp. 886--905, 1996.
// There appears to be a plethora of publicly available preprints.
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
  const Int num_orig_vert = quotient_graph.num_original_vertices;

  // Initialize the shrinkable cached binary tree of external degrees.
  ShrinkableCachedBinaryTree<Int> external_degrees;
  {
    std::vector<Int> external_degrees_vec(num_orig_vert);
    for (Int source = 0; source < num_orig_vert; ++source) {
      external_degrees_vec[source] =
          quotient_graph.adjacency_lists[source].size();
    }
    external_degrees.Initialize(external_degrees_vec);
  }

  // Buffers used to store a temporary integer vector.
  // They are declared here to avoid unnecessary garbage collection.
  std::vector<Int> temp_int_vec0, temp_int_vec1;

  // Eliminate the variables.
  MinimumDegreeAnalysis analysis(num_orig_vert);
  while (quotient_graph.num_eliminated_vertices < num_orig_vert) {
    // Retrieve a variable with minimal (approximate) external degree.
    const std::pair<Int, Int> pivot_pair = external_degrees.MinimalEntry(); 
    const Int pivot = pivot_pair.first;

    // Compute the structure of the pivot:
    //   L_p := (A_p \cup (\cup_{e in E_p} L_e)) \ p.
    //
    // TODO(Jack Poulson): Experiment with several different union-finding
    // algorithms. It isn't clear that repeated calls to std::set_union has
    // optimal time complexity.
    FilterSet(
      quotient_graph.adjacency_lists[pivot],
      quotient_graph.supervariables[pivot],
      &quotient_graph.structures[pivot]);
    for (const Int& element : quotient_graph.element_lists[pivot]) {
      FilterSet(
        quotient_graph.structures[element],
        quotient_graph.supervariables[pivot],
        &temp_int_vec0);

      temp_int_vec1 = quotient_graph.structures[pivot];
      MergeSets(
        temp_int_vec1,
        temp_int_vec0,
        &quotient_graph.structures[pivot]);
    }

    // Update the adjacency lists.
    for (const Int& i : quotient_graph.structures[pivot]) {
      // Remove redundant adjacency entries:
      //   A_i := (A_i \ L_p) \ supervar(p).
      FilterSet(
          quotient_graph.adjacency_lists[i],
          quotient_graph.structures[pivot],
          &temp_int_vec0);
      FilterSet(
          temp_int_vec0,
          quotient_graph.supervariables[pivot],
          &quotient_graph.adjacency_lists[i]);
    }

    // Update the element lists.
    for (const Int& i : quotient_graph.structures[pivot]) {
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
    std::unordered_map<Int, Int> external_structure_cardinalities;
    if (degree_type == kAmestoyExternalDegree) {
      external_structure_cardinalities = ExternalStructureCardinalities(
          quotient_graph, pivot);
    }

    // Compute the external degree approximations of the supervariables
    // adjacent to the current pivot.
    for (const Int& i : quotient_graph.structures[pivot]) {
      // Compute the external degree (or approximation) of supervariable i:
      //   d_i := |A_i \ supervar(i)| + |(\cup_{e in E_i} L_e) \ supervar(i)|.
      const Int old_external_degree = external_degrees.GetValue(i); 
      const Int external_degree = ExternalDegree(
        quotient_graph, i, pivot, old_external_degree,
        external_structure_cardinalities, degree_type);
      external_degrees.SetValue(i, external_degree);
    }

    // Fill a set of buckets for the hashes of the supervariables adjacent to
    // the current pivot.
    const Int struct_size = quotient_graph.structures[pivot].size();
    std::vector<Int> bucket_list(struct_size);
    std::unordered_map<Int, std::unordered_set<Int>> variable_hash_map;
    for (Int i_index = 0; i_index < struct_size; ++i_index) {
      const Int i = quotient_graph.structures[pivot][i_index];

      // Append this principal variable to its hash bucket.
      const Int bucket = AshcraftVariableHash(
        quotient_graph.num_orig_vertices,
        quotient_graph.adjacency_list[i],
        quotient_graph.element_list[i]);
      bucket_list[i_index] = bucket;
      variable_hash_map[bucket].push_back(i_index);
    }

    // Supervariable detection.
    for (Int i_index = 0; i_index < pivot_struct_size; ++i_index) {
      const Int i = quotient_graph.structures[pivot][i_index];
      const Int bucket = bucket_list[i_index];
      for (const Int& j_index : variable_hash_map[bucket]) {
        if (j_index <= i_index) {
          continue;
        }
        const Int j = quotient_graph.structures[pivot][j_index];
        const bool indistinguishable =
            SupernodesAreQuotientIndistinguishable(quotient_graph, i, j);
        if (indistinguishable) {
          // Absorb supervar(j) into supervar(i). 
          {
            temp_int_vec0 = quotient_graph.supervariables[i];
            MergeSets(
                temp_int_vec0,
                quotient_graph.supervariables[j],
                &quotient_graph.supervariables[i]);
          }
          external_degrees.UpdateValue(
              i, -quotient_graph.supervariables[j].size());
          external_degrees.DisableIndex(j);
#ifdef QUOTIENT_DEBUG
          quotient_graph.variables.erase(j);
#endif
          SwapClearVector(quotient_graph.supervariables[j]);
          SwapClearVector(quotient_graph.adjacency_lists[j]);
          SwapClearVector(quotient_graph.element_lists[j]);
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
    SwapClearVector(quotient_graph.adjacency_lists[pivot]);
    SwapClearVector(quotient_graph.element_lists[pivot]);
    quotient_graph.elimination_order.push_back(pivot);
    quotient_graph.num_eliminated_vertices +=
        quotient_graph.supervariables[pivot].size();
  }

  analysis.supervariables = quotient_graph.supervariables;
  analysis.structures = quotient_graph.structures;
  return analysis;
}

} // namespace quotient

#endif // ifndef QUOTIENT_MINIMUM_DEGREE_H_
