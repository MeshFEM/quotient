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
inline std::size_t AshcraftHash(
    Int num_vertices,
    const std::vector<Int>& adjacency_list,
    const std::vector<Int>& element_list) {
  std::size_t result = 0;
  for (const Int& index : adjacency_list) {
    result = (result + index) % (num_vertices - 1);
  }
  for (const Int& index : element_list) {
    result = (result + index) % (num_vertices - 1);
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
    Int i,
    Int j,
    const std::vector<Int>& adjacency_list_i,
    const std::vector<Int>& adjacency_list_j,
    const std::vector<Int>& element_list_i,
    const std::vector<Int>& element_list_j) {
  const Int adj_i_size = adjacency_list_i.size();
  const Int adj_j_size = adjacency_list_j.size();
  const Int elem_i_size = element_list_i.size();
  const Int elem_j_size = element_list_j.size();

  // Early exit if the set cardinalities disagree.
  if (adj_i_size + elem_i_size != adj_j_size + elem_j_size) {
    return false;
  }

  // Store indices into the adjacency list, element list, and singleton.
  Int[3] i_heads{0, 0, 0};
  Int[3] j_heads{0, 0, 0};
  for (Int index = 0; index < adj_i_size + elem_i_size + 1; ++index) {
    // Test for equality with adjacency_list_i[i_heads[0]].
    if (i_heads[0] < adj_i_size) {
      const Int i_candidate = adjacency_list_i[i_heads[0]];
      if (j_heads[0] < adj_j_size &&
          i_candidate == adjacency_list_j[j_heads[0]]) {
        ++i_heads[0];
        ++j_heads[0];
        continue;
      }
      if (j_heads[1] < elem_j_size &&
                 i_candidate == element_list_j[j_heads[1]]) {
        ++i_heads[0];
        ++j_heads[1];
        continue;
      }
      if (j_heads[2] == 0 && i_candidate == j) {
        ++i_heads[0];
        ++j_heads[2];
        continue;
      }
    }

    // Test for equality with element_list_i[i_heads[1]].
    if (i_heads[1] < elem_i_size) {
      const Int i_candidate = element_list_i[i_heads[1]];
      if (j_heads[0] < adj_j_size &&
          i_candidate == adjacency_list_j[j_heads[0]]) {
        ++i_heads[1];
        ++j_heads[0];
        continue;
      }
      if (j_heads[1] < elem_j_size &&
                 i_candidate == element_list_j[j_heads[1]]) {
        ++i_heads[1];
        ++j_heads[1];
        continue;
      }
      if (j_heads[2] == 0 && i_candidate == j) {
        ++i_heads[1];
        ++j_heads[2];
        continue;
      }
    }

    // Test for equality with i.
    if (i_heads[2] == 0) {
      const Int i_candidate = i;
      if (j_heads[0] < adj_j_size &&
          i_candidate == adjacency_list_j[j_heads[0]]) {
        ++i_heads[2];
        ++j_heads[0];
        continue;
      }
      if (j_heads[1] < elem_j_size &&
                 i_candidate == element_list_j[j_heads[1]]) {
        ++i_heads[2];
        ++j_heads[1];
        continue;
      }
      if (j_heads[2] == 0 && i_candidate == j) {
        ++i_heads[2];
        ++j_heads[2];
        continue;
      }
    }

    return false;
  }

#ifdef QUOTIENT_DEBUG
  if (i_heads[0] != adjacency_list_i.size() ||
      i_heads[1] != element_list_i.size() ||
      i_heads[2] != 1) {
    std::cerr << "i heads did not match list sizes." << std::endl;
  }
  if (j_heads[0] != adjacency_list_j.size() ||
      j_heads[1] != element_list_j.size() ||
      j_heads[2] != 1) {
    std::cerr << "j heads did not match list sizes." << std::endl;
  }
#endif

  return true;
}

inline Int ExactExternalDegree(
    const std::vector<Int>& supervariable,
    const std::vector<Int>& adjacency_list,
    const std::vector<Int>& element_list,
    const std::vector<std::vector<Int>>& structures) {
  // TODO(Jack Poulson)
}

inline Int AmestoyExternalDegree(
    const std::vector<Int>& supervariable,
    const std::vector<Int>& adjacency_list,
    const std::vector<Int>& element_list,
    const std::vector<std::vector<Int>>& structures,
    Int num_vertices_left,
    Int old_external_degree) {
  // TODO(Jack Poulson)
}

// Returns the external degree approximation of Gilbert, Moler, and Schreiber,
//   \hat{d_i} = |A_i \ supervar(i)|  + \sum_{e in E_i} |L_e \ supervar(i)|.
inline Int GilbertExternalDegree(
    const std::vector<Int>& supervariable,
    const std::vector<Int>& adjacency_list,
    const std::vector<Int>& element_list,
    const std::vector<std::vector<Int>>& structures) {
  Int external_degree = adjacency_list.size() -
      CardinalityOfIntersection(adjacency_list, supervariable);
  for (const Int& element : element_list) {
    external_degree += structures[element] -
        CardinalityOfIntersection(structures[element], supervariable);
  }
  return external_degree;
}

// Returns the external degree approximation of Ashcraft, Eisenstat, and Lucas:
//   \tilde{d_i} = d_i if |E_i| = 2, \hat{d_i} otherwise.
inline Int AshcraftExternalDegree(
    const std::vector<Int>& supervariable,
    const std::vector<Int>& adjacency_list,
    const std::vector<Int>& element_list,
    const std::vector<std::vector<Int>>& structures) {
  if (element_list.size() == 2) {

  }
  return GilbertExternalDegree(
      supervariable, adjacency_list, element_list, structures);
}

// Returns (an approximation of) the external degree of a given supervariable.
inline Int ExternalDegree(
    const std::vector<Int>& supervariable,
    const std::vector<Int>& adjacency_list,
    const std::vector<Int>& element_list,
    const std::vector<std::vector<Int>>& structures,
    Int num_vertices_left,
    Int old_external_degree,
    ExternalDegreeType degree_type) {
  switch(degree_type) {
    case kExactExternalDegree: { 
      return ExactExternalDegree(
          supervariable, adjacency_list, element_list, structures);
    }
    case kAmestoyExternalDegree: {
      return AmestoyExternalDegree(
          supervariable, adjacency_list, element_list, structures,
          num_vertices_left, old_external_degree);
    }
    case kAshcraftExternalDegree: { 
      return AshcraftExternalDegree(
          supervariable, adjacency_list, element_list, structures);
    }
    case kGilbertExternalDegree: {
      return GilbertExternalDegree(
          supervariable, adjacency_list, element_list, structures);
    }
  }
}

struct MinimumDegreeAnalysis {
  std::vector<std::vector<Int>> supervariables;
  std::vector<Int> elimination_order;
  std::vector<std::vector<Int>> structures;

  Analysis(Int num_vertices) {
    supervariables.resize(num_vertices);
    for (Int source = 0; source < num_vertices; ++source) {
      supervariables[source].push_back(source);
    }
    elimination_order.reserve(num_vertices);
    structures.resize(num_vertices);
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
  const Int num_vertices = graph.NumSources();
  const std::vector<GraphEdge>& graph.Edges();

#ifdef QUOTIENT_DEBUG
  // Initialize the (principal) variable and element lists.
  std::set<Int> variables;
  for (Int source = 0; source < num_vertices; ++source) {
    variables.insert(source);
  }
  std::set<Int> elements;
#endif

  // Initialize the adjacenty and lists from the original graph.
  std::vector<std::vector<Int>> adjacency_lists(num_vertices);
  std::vector<std::vector<Int>> element_lists(num_vertices);
  for (Int source = 0; source < num_vertices; ++source) {
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

  // Initialize the shrinkable cached binary tree of external degrees.
  ShrinkableCachedBinaryTree<Int> external_degrees;
  {
    std::vector<Int> external_degrees_vec(num_vertices);
    for (Int source = 0; source < num_vertices; ++source) {
      external_degrees_vec[source] = adjacency_lists[source].size();
    }
    external_degrees.Initialize(external_degrees_vec);
  }

  // Buffers used to store a temporary integer vector.
  // They are declared here to avoid unnecessary garbage collection.
  std::vector<Int> temp_int_vec0, temp_int_vec1;

  // Eliminate the variables.
  MinimumDegreeAnalysis analysis(num_vertices);
  Int num_eliminated = 0;
  while (num_eliminated < num_vertices) {
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
      adjacency_lists[pivot],
      analysis.supervariables[pivot],
      &analysis.structures[pivot]);
    for (const Int& element : element_lists[pivot]) {
      FilterSet(
        analysis.structures[element],
        analysis.supervariables[pivot],
        &temp_int_vec0);

      temp_int_vec1 = analysis.structures[pivot];
      MergeSets(
        temp_int_vec1,
        temp_int_vec0,
        &analysis.structures[pivot]);
    }

    // A set of buckets for hashes of the supervariables adjacent to the
    // current pivot.
    std::vector<Int> bucket_list;
    std::unordered_map<Int, std::unordered_set<Int>> variable_hash_map;

    // Update the adjacency lists, element lists, and external degrees.
    const Int struct_size = analysis.structures[pivot].size();
    for (Int i_index = 0; i_index < struct_size; ++i_index) {
      const Int i = analysis.structures[pivot][i_index];

      // Remove redundant adjacency entries:
      //   A_i := (A_i \ L_p) \ supervar(p).
      FilterSet(
          adjacency_lists[i],
          analysis.structures[pivot],
          &temp_int_vec0);
      FilterSet(
          temp_int_vec0,
          analysis.supervariables[pivot],
          &adjacency_lists[i]);

      // Element absorption:
      //   E_i := (E_i \ E_p) \cup {p}.
      temp_int_vec0 = element_lists[i];
      FilterSet(
          temp_int_vec0,
          element_lists[p],
          &element_lists[i]);
      InsertEntryIntoSet(pivot, &element_lists[i]);

      // Compute external degree (or approximation of):
      //   d_i := |A_i \ supervar(i)| + |(\cup_{e in E_i} L_e) \ supervar(i)|.
      const Int old_external_degree = external_degrees.GetValue(i); 
      const Int num_vertices_left = num_vertices - num_eliminated;
      const Int external_degree = ExternalDegree(
        analysis.supervariables[i],
        adjacency_lists[i],
        element_lists[i],
        analysis.structures,
        num_vertices_left,
        old_external_degree,
        degree_type);
      external_degrees.SetValue(i, external_degree);

      // Append this principal variable to its hash bucket.
      const Int bucket = AshcraftHash(
        num_variables, adjacency_list[i], element_list[i]);
      bucket_list[i_index] = bucket;
      variable_hash_map[bucket].push_back(i_index);
    }

    // Supervariable detection.
    for (Int i_index = 0; i_index < pivot_struct_size; ++i_index) {
      const Int i = analysis.structures[pivot][i_index];
      const Int bucket = bucket_list[i_index];
      for (const Int& j_index : variable_hash_map[bucket]) {
        if (j_index <= i_index) {
          continue;
        }
        const Int j = analysis.structures[pivot][j_index];
        const bool indistinguishable =
            SupernodesAreQuotientIndistinguishable(
                i, j,
                adjacency_lists[i], adjacency_lists[j],
                element_lists[i], element_lists[j]);
        if (indistinguishable) {
          // Absorb supervar(j) into supervar(i). 
          {
            temp_int_vec0 = analysis.supervariables[i];
            MergeSets(
                temp_int_vec0,
                analysis.supervariables[j],
                &analysis.supervariables[i]);
          }
          external_degrees.UpdateValue(i, -analysis.supervariables[j].size());
          external_degrees.DisableIndex(j);
#ifdef QUOTIENT_DEBUG
          variables.erase(j);
#endif
          SwapClearVector(analysis.supervariables[j]);
          SwapClearVector(adjacency_lists[j]);
          SwapClearVector(element_lists[j]);
        }
      }
    }

    // Convert pivot from a variable to an element.
#ifdef QUOTIENT_DEBUG
    // Update the element list:
    //   \bar{V} := (\bar{V} \cup {p}) \ E_p
    elements.insert(pivot);
    for (const Int& element : element_lists[pivot]) {
      elements.erase(element);
    }
    // Update the variable list:
    //   V := V \ {p}.
    variables.erase(pivot);
#endif
    SwapClearVector(adjacency_lists[pivot]);
    SwapClearVector(element_lists[pivot]);
    analysis.elimination_order.push_back(pivot);
    num_eliminated += analysis.supervariables[pivot].size();
  }

  return analysis;
}

} // namespace quotient

#endif // ifndef QUOTIENT_MINIMUM_DEGREE_H_
