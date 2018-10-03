/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_EXTERNAL_DEGREE_H_
#define QUOTIENT_EXTERNAL_DEGREE_H_

#include <iostream>
#include <unordered_map>
#include <vector>

#include "quotient/config.hpp"
#include "quotient/quotient_graph.hpp"
#include "quotient/set_utils.hpp"

namespace quotient {

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

} // namespace quotient

#endif // ifndef QUOTIENT_EXTERNAL_DEGREE_H_
