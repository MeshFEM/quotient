/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_EXTERNAL_DEGREE_IMPL_H_
#define QUOTIENT_EXTERNAL_DEGREE_IMPL_H_

#include <iostream>
#include <vector>

#include "quotient/config.hpp"
#include "quotient/quotient_graph.hpp"
#include "quotient/set_utils.hpp"
#include "quotient/external_degree.hpp"

namespace quotient {

namespace external_degree {

namespace exact {

// Computes the exact external degree of supernode i using Eq. (2) of [ADD-96].
//   d_i = |A_i \ supernode(i)| + |(\cup_{e in E_i) L_e) \ supernode(i)|.
//
// This version is specialized to the case where |E_i| = 2.
inline Int TwoElements(const QuotientGraph& graph, Int i) {
#ifdef QUOTIENT_DEBUG
  if (graph.element_lists[i].size() != 2) {
    std::cerr << "Expected exactly two elements in list." << std::endl;
  }
#endif
  // Add the cardinality of A_i \ supernode(i).
  Int degree = SizeOfDifference(graph.adjacency_lists[i], graph.supernodes[i]);

  // Add the cardinality of (\cup_{e in E_i} L_e) \ supernode(i).
  const Int element0 = graph.element_lists[i][0];
  const Int element1 = graph.element_lists[i][1];
  degree += SizeOfBlacklistedUnion(
    graph.structures[element0],
    graph.structures[element1],
    graph.supernodes[i]);

  return degree;
}

// Computes the exact external degree of supernode i using Eq. (2) of [ADD-96].
//   d_i = |A_i \ supernode(i)| + |(\cup_{e in E_i) L_e) \ supernode(i)|.
//
// This version is specialized to the case where |E_i| is small.
inline Int FewElements(const QuotientGraph& graph, Int i) {
  // Add the cardinality of A_i \ supernode(i).
  Int degree = SizeOfDifference(graph.adjacency_lists[i], graph.supernodes[i]);

  // Add the cardinality of (\cup_{e in E_i} L_e) \ supernode(i).
  std::vector<Int> filtered_struct;
  std::vector<Int> element_struct_union;
  std::vector<Int> temp_int_vec;
  for (const Int& element : graph.element_lists[i]) {
    FilterSet(graph.structures[element], graph.supernodes[i],
              &filtered_struct);
    temp_int_vec = element_struct_union;
    MergeSets(temp_int_vec, filtered_struct, &element_struct_union);
  }
  degree += element_struct_union.size();

  return degree;
}

// Computes the exact external degree of supernode i using Eq. (2) of [ADD-96].
//   d_i = |A_i \ supernode(i)| + |(\cup_{e in E_i) L_e) \ supernode(i)|.
//
// This version is specialized to the case where |E_i| is not small.
inline Int ManyElements(const QuotientGraph& graph, Int i) {
  // Add the cardinality of A_i \ supernode(i).
  Int degree = SizeOfDifference(graph.adjacency_lists[i], graph.supernodes[i]);

  // Count the number of (non-unique) filtered entries and compute the scan
  // of the offsets.
  Int num_non_unique = 0;
  std::vector<Int> offsets;
  offsets.reserve(graph.element_lists[i].size());
  for (const Int& element : graph.element_lists[i]) {
    offsets.push_back(num_non_unique);
    num_non_unique += SizeOfDifference(
        graph.structures[element], graph.supernodes[i]);
  }

  // Fill the unsorted, non-unique list of elements.
  //
  // TODO(Jack Poulson): Use OpenMP to parallelize this loop.
  std::vector<Int> union_of_structures(num_non_unique);
  for (std::size_t index = 0; index < graph.element_lists[i].size(); ++index) {
    const Int& element = graph.element_lists[i][index];
    const std::vector<Int>& structure = graph.structures[element];
    std::set_difference(
      structure.begin(), structure.end(),
      graph.supernodes[i].begin(), graph.supernodes[i].end(),
      union_of_structures.begin() + offsets[index]);
  }

  // Sort the non-unique list of elements.
  std::sort(union_of_structures.begin(), union_of_structures.end());

  // Add the number of unique entries onto the external degree.
  //
  // TODO(Jack Poulson): Use OpenMP to parallelize this.
  Int num_unique = 0;
  Int last_entry = -1;
  for (const Int& entry : union_of_structures) {
    if (entry != last_entry) {
      ++num_unique;
      last_entry = entry;
    }
  }
  degree += num_unique;

  return degree;
}

}  // namespace exact

// Computes the exact external degree of supernode i using Eq. (2) of [ADD-96].
//   d_i = |A_i \ supernode(i)| + |(\cup_{e in E_i) L_e) \ supernode(i)|.
inline Int Exact(const QuotientGraph& graph, Int i) {
  const Int kElementThreshold = 4;

  if (graph.element_lists[i].size() == 2) {
    return exact::TwoElements(graph, i);
  }
  if (graph.element_lists[i].size() < kElementThreshold) {
    return exact::FewElements(graph, i);
  }
  return exact::ManyElements(graph, i);
}

// Computes an approximation of the external degree of supernode i using Eq. (4)
// of [ADD-96].
inline Int Amestoy(
    const QuotientGraph& graph,
    Int i,
    Int pivot,
    const std::vector<Int>& external_structure_sizes) {
  const Int num_vertices_left =
      graph.num_original_vertices - graph.num_eliminated_vertices;
  const Int old_degree = graph.external_degree_heap.Value(i); 

  // Note that this usage of 'external' refers to |L_p \ supernode(i)| and not
  // |L_e \ L_p|, as is the case for 'external_structure_sizes'.
  const Int external_pivot_structure_size = SizeOfDifference(
      graph.structures[pivot], graph.supernodes[i]);

  const Int bound0 = num_vertices_left;
  const Int bound1 = old_degree + external_pivot_structure_size;

  // bound_2 = |A_i \ supernode(i)| + |L_p \ supernode(i)| + 
  //           \sum_{e in E_i \ {p}} |L_e \ L_p|.
  Int bound2 =
    SizeOfDifference(graph.adjacency_lists[i], graph.supernodes[i]) +
    SizeOfDifference(graph.structures[pivot], graph.supernodes[i]);
  for (const Int& element : graph.element_lists[i]) {
    if (element == pivot) {
      continue;
    }
    if (external_structure_sizes[element] >= 0) {
      bound2 += external_structure_sizes[element];
    } else {
      bound2 += graph.structures[element].size();
    }
  }

  return std::min(bound0, std::min(bound1, bound2));
}

// Returns the external degree approximation of Gilbert, Moler, and Schreiber,
//   \hat{d_i} = |A_i \ supernode(i)|  + \sum_{e in E_i} |L_e \ supernode(i)|.
inline Int Gilbert(const QuotientGraph& graph, Int i) {
  Int degree = SizeOfDifference(graph.adjacency_lists[i], graph.supernodes[i]);
  for (const Int& element : graph.element_lists[i]) {
    degree += SizeOfDifference(graph.structures[element], graph.supernodes[i]);
  }
  return degree;
}

// Returns the external degree approximation of Ashcraft, Eisenstat, and Lucas:
//   \tilde{d_i} = d_i if |E_i| = 2, \hat{d_i} otherwise.
inline Int Ashcraft(const QuotientGraph& graph, Int i) {
  if (graph.element_lists[i].size() == 2) {
    return Exact(graph, i);
  }
  return Gilbert(graph, i);
}

}  // namespace external_degree

inline Int ExternalDegree(
    const QuotientGraph& graph,
    Int i,
    Int pivot,
    const std::vector<Int>& external_structure_sizes,
    ExternalDegreeType degree_type) {
  Int degree = -1;
  switch(degree_type) {
    case kExactExternalDegree: { 
      degree = external_degree::Exact(graph, i);
      break;
    }
    case kAmestoyExternalDegree: {
      degree = external_degree::Amestoy(
          graph, i, pivot, external_structure_sizes);
      break;
    }
    case kAshcraftExternalDegree: { 
      degree = external_degree::Ashcraft(graph, i);
      break;
    }
    case kGilbertExternalDegree: {
      degree = external_degree::Gilbert(graph, i);
      break;
    }
  }
  return degree;
}

} // namespace quotient

#endif // ifndef QUOTIENT_EXTERNAL_DEGREE_IMPL_H_
