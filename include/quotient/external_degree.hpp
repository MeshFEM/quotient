/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_EXTERNAL_DEGREE_H_
#define QUOTIENT_EXTERNAL_DEGREE_H_

#include <vector>

#include "quotient/config.hpp"
#include "quotient/quotient_graph.hpp"

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
  //         \sum_{e in E_i \ {p}} |L_e \ L_p|).
  kAmestoyExternalDegree,

  // In the notation of [ADD-96],
  //   \tilde{d_i} = d_i if |E_i| = 2, \hat{d_i}, otherwise.
  kAshcraftExternalDegree,

  // In the notation of [ADD-96],
  //   \hat{d_i} = |A_i \ supernode(i)| + \sum_{e in E_i} |L_e \ supernode(i)|.
  kGilbertExternalDegree,
};


// Returns (an approximation of) the external degree of a given supervariable.
Int ExternalDegree(
    const QuotientGraph& graph,
    Int i,
    Int pivot,
    const std::vector<int>& pivot_structure_mask,
    const std::vector<Int>& external_structure_sizes,
    ExternalDegreeType degree_type);

} // namespace quotient

#include "quotient/external_degree-impl.hpp"

#endif // ifndef QUOTIENT_EXTERNAL_DEGREE_H_
