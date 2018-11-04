/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_DEGREE_TYPE_H_
#define QUOTIENT_DEGREE_TYPE_H_

namespace quotient {

// A sequence of external degree approximations of decreasing accuracy.
// Please see Theorem 4.1 of [ADD-96] for accuracy guarantees.
enum DegreeType {
  // In the notation of [ADD-96],
  //   d_i = |A_i \ supernode(i)| + |(\cup_{e in E_i} L_e) \ supernode(i)|.
  kExactDegree,

  // In the notation of [ADD-96],
  //   \bar{d_i}^k = min(
  //     n - k,
  //     \bar{d_i}^{k - 1} + |L_p \ supernode(i)|, 
  //     |A_i \ supernode(i)| + |L_p \ supernode(i)| +
  //         \sum_{e in E_i \ {p}} |L_e \ L_p|).
  kAmestoyDegree,

  // In the notation of [ADD-96],
  //   \tilde{d_i} = d_i if |E_i| = 2, \hat{d_i}, otherwise.
  kAshcraftDegree,

  // In the notation of [ADD-96],
  //   \hat{d_i} = |A_i \ supernode(i)| + \sum_{e in E_i} |L_e \ supernode(i)|.
  kGilbertDegree,
};

} // namespace quotient

#endif // ifndef QUOTIENT_DEGREE_TYPE_H_
