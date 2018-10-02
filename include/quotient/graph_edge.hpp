/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_GRAPH_EDGE_H_
#define QUOTIENT_GRAPH_EDGE_H_

#include "quotient/config.hpp"

namespace quotient {

// A simple structure for storing graph edges.
struct GraphEdge {
  // The index of the source node where the edge begins.
  Int source;

  // The index of the target node where the edge ends.
  Int target;

  // The default constructor.
  GraphEdge() { }

  // The fully-specified constructor.
  GraphEdge(Int source_val, Int target_val)
  : source(source_val), target(target_val) { }

  // A lexicographic less-than operator for GraphEdge that can be passed into
  // STL routines as a comparison function.
  static bool LessThan(const GraphEdge& lhs, const GraphEdge& rhs) {
    return (lhs.source < rhs.source) ||
        (lhs.source == rhs.source && lhs.target < rhs.target);
  }

  // A lexicographic greater-than operator for GraphEdge that can be passed
  // into STL routines as a comparison function.
  static bool GreaterThan(const GraphEdge& lhs, const GraphEdge& rhs) {
    return (lhs.source > rhs.source) ||
        (lhs.source == rhs.source && lhs.target > rhs.target);
  }
};

// The lexicographic less-than comparison operator for GraphEdge.
inline bool operator<(const GraphEdge& lhs, const GraphEdge& rhs) {
  return GraphEdge::LessThan(lhs, rhs);
}

// The lexicographic less-than-or-equal comparison operator for GraphEdge.
inline bool operator<=(const GraphEdge& lhs, const GraphEdge& rhs) {
  return !GraphEdge::GreaterThan(lhs, rhs);
}

// The lexicographic greater-than comparison operator for GraphEdge.
inline bool operator>(const GraphEdge& lhs, const GraphEdge& rhs) {
  return GraphEdge::GreaterThan(lhs, rhs);
}

// The lexicographic greater-than-or-equal comparison operator for GraphEdge.
inline bool operator>=(const GraphEdge& lhs, const GraphEdge& rhs) {
  return !GraphEdge::LessThan(lhs, rhs);
}

// An equality check for two GraphEdge instances.
inline bool operator==(const GraphEdge& lhs, const GraphEdge& rhs) {
  return lhs.source == rhs.source && lhs.target == rhs.target;
}

// An inequality check for two GraphEdge instances.
inline bool operator!=(const GraphEdge& lhs, const GraphEdge& rhs) {
  return !(lhs == rhs);
}

} // namespace quotient

#endif // ifndef QUOTIENT_GRAPH_EDGE_H_
