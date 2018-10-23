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

// Computes the exact external degree of supernode i using Eq. (2) of [ADD-96].
//   d_i = |A_i \ supernode(i)| + |(\cup_{e in E_i) L_e) \ supernode(i)|.
inline Int Exact(const QuotientGraph& graph, Int i, std::vector<int>* mask) {
  // Add the cardinality of A_i \ supernode(i).
  Int degree = 0;
  for (const Int& index : graph.adjacency_lists[i]) {
    degree += graph.supernode_sizes[index];
  }

  // Add on the number of unique entries in the structures of the element lists
  // that are outside supernode(i).
  for (const Int& element : graph.element_lists[i]) {
    for (const Int& j : graph.structures[element]) {
      const Int head = graph.head_index[j];
      if ((*mask)[j] || head == i || graph.supernode_sizes[head] < 0) {
        continue;
      }
      ++degree;
      (*mask)[j] = 1;
    }
  }

  // Clear the mask.
  // NOTE: We allow (*mask)[i] to be set to zero to avoid branching.
  for (const Int& element : graph.element_lists[i]) {
    for (const Int& j : graph.structures[element]) {
      (*mask)[j] = 0;
    }
  }

  return degree;
}

// Computes an approximation of the external degree of supernode i using Eq. (4)
// of [ADD-96].
//
// This routine is made slightly more accurate by subtracting the size of
// supernode(i) from bound0.
inline Int Amestoy(
    const QuotientGraph& graph,
    Int i,
    Int pivot,
    const std::vector<int>& pivot_structure_mask,
    const std::vector<Int>& external_structure_sizes) {
  const Int num_vertices_left =
      graph.num_original_vertices - graph.num_eliminated_vertices;
  const Int old_degree = graph.degree_lists.degrees[i];

  // Note that this usage of 'external' refers to |L_p \ supernode(i)| and not
  // |L_e \ L_p|, as is the case for 'external_structure_sizes'.
  Int external_pivot_structure_size = graph.structures[pivot].size();
  if (pivot_structure_mask[i]) {
    external_pivot_structure_size -= graph.supernode_sizes[i];
  }

  const Int bound0 = num_vertices_left - graph.supernode_sizes[i];
  const Int bound1 = old_degree + external_pivot_structure_size;

  // bound_2 = |A_i \ supernode(i)| + |L_p \ supernode(i)| +
  //           \sum_{e in E_i \ {p}} |L_e \ L_p|.
  Int bound2 = 0;
  for (const Int& index : graph.adjacency_lists[i]) {
    bound2 += graph.supernode_sizes[index];
  }
  bound2 += external_pivot_structure_size;
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
//
// We slightly modify this formula to ensure that the estimated degree is
// at most
//   (num_original_vertices - num_eliminated_vertices) - size(supernode(i)).
inline Int Gilbert(const QuotientGraph& graph, Int i) {
  Int degree = 0;
  for (const Int& index : graph.adjacency_lists[i]) {
    degree += graph.supernode_sizes[index];
  }
  for (const Int& element : graph.element_lists[i]) {
    for (const Int& j : graph.structures[element]) {
      const Int head = graph.head_index[j];
      if (head == i || graph.supernode_sizes[head] < 0) {
        continue;
      }
      ++degree;
    }
  }
  const Int num_vertices_left =
      graph.num_original_vertices - graph.num_eliminated_vertices;
  return std::min(degree, num_vertices_left - graph.supernode_sizes[i]);
}

// Returns the external degree approximation of Ashcraft, Eisenstat, and Lucas:
//   \tilde{d_i} = d_i if |E_i| = 2, \hat{d_i} otherwise.
inline Int Ashcraft(
    const QuotientGraph& graph, Int i, std::vector<int>* exact_degree_mask) {
  if (graph.element_lists[i].size() == 2) {
    return Exact(graph, i, exact_degree_mask);
  }
  return Gilbert(graph, i);
}

}  // namespace external_degree

inline Int ExternalDegree(
    const QuotientGraph& graph,
    Int i,
    Int pivot,
    const std::vector<int>& pivot_structure_mask,
    const std::vector<Int>& external_structure_sizes,
    ExternalDegreeType degree_type,
    std::vector<int>* exact_degree_mask) {
  Int degree = -1;
  switch(degree_type) {
    case kExactExternalDegree: {
      degree = external_degree::Exact(graph, i, exact_degree_mask);
      break;
    }
    case kAmestoyExternalDegree: {
      degree = external_degree::Amestoy(
          graph, i, pivot, pivot_structure_mask, external_structure_sizes);
      break;
    }
    case kAshcraftExternalDegree: {
      degree = external_degree::Ashcraft(graph, i, exact_degree_mask);
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
