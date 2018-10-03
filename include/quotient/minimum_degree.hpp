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
#include <vector>

#include "quotient/config.hpp"
#include "quotient/coordinate_graph.hpp"
#include "quotient/external_degree.hpp"
#include "quotient/quotient_graph.hpp"
#include "quotient/random_access_heap.hpp"

namespace quotient {

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
