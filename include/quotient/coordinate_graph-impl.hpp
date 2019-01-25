/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_COORDINATE_GRAPH_IMPL_H_
#define QUOTIENT_COORDINATE_GRAPH_IMPL_H_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

#include "quotient/coordinate_graph.hpp"

namespace quotient {

template <typename T>
void SwapClearVector(std::vector<T>* vec) {
  std::vector<T>().swap(*vec);
}

template <typename T>
void EraseDuplicatesInSortedVector(std::vector<T>* vec) {
  vec->erase(std::unique(vec->begin(), vec->end()), vec->end());
}

template <typename T>
void EraseDuplicatesInSortedVector(Buffer<T>* vec) {
  T* new_end = std::unique(vec->begin(), vec->end());
  const Int new_size = std::distance(vec->begin(), new_end);
  vec->Resize(new_size);
}

inline CoordinateGraph::CoordinateGraph() : num_sources_(0), num_targets_(0) {}

inline CoordinateGraph::CoordinateGraph(const CoordinateGraph& graph) {
  if (&graph == this) {
    return;
  }
  *this = graph;
}

inline CoordinateGraph& CoordinateGraph::operator=(
    const CoordinateGraph& graph) {
  if (&graph != this) {
    num_sources_ = graph.num_sources_;
    num_targets_ = graph.num_targets_;
    edges_ = graph.edges_;
    source_edge_offsets_ = graph.source_edge_offsets_;

    edges_to_add_ = graph.edges_to_add_;
    edges_to_remove_ = graph.edges_to_remove_;
  }
  return *this;
}

inline std::unique_ptr<CoordinateGraph> CoordinateGraph::FromMatrixMarket(
    const std::string& filename, bool skip_explicit_zeros, EntryMask mask) {
  std::unique_ptr<CoordinateGraph> result;
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Could not open " << filename << std::endl;
    return result;
  }

  // Fill the description of the Matrix Market data.
  MatrixMarketDescription description;
  if (!ReadMatrixMarketDescription(file, &description)) {
    return result;
  }

  result.reset(new CoordinateGraph);
  if (description.format == kMatrixMarketFormatArray) {
    // Read the size of the matrix.
    Int num_rows, num_columns;
    if (!ReadMatrixMarketArrayMetadata(description, file, &num_rows,
                                       &num_columns)) {
      result.reset();
      return result;
    }

    // Fill a fully-connected graph.
    result->AsymmetricResize(num_rows, num_columns);
    result->ReserveEdgeAdditions(num_rows * num_columns);
    for (Int row = 0; row < num_rows; ++row) {
      for (Int column = 0; column < num_columns; ++column) {
        result->QueueEdgeAddition(row, column);
      }
    }
    result->FlushEdgeQueues();
    return result;
  }

  // Read in the number of matrix dimensions and the number of entries specified
  // in the file.
  Int num_rows, num_columns, num_entries;
  if (!ReadMatrixMarketCoordinateMetadata(description, file, &num_rows,
                                          &num_columns, &num_entries)) {
    result.reset();
    return result;
  }

  // Fill in the edges.
  Int num_skipped_entries = 0;
  const Int num_entries_bound =
      description.symmetry == kMatrixMarketSymmetryGeneral ? num_entries
                                                           : 2 * num_entries;
  result->AsymmetricResize(num_rows, num_columns);
  result->ReserveEdgeAdditions(num_entries_bound);
  for (Int edge_index = 0; edge_index < num_entries; ++edge_index) {
    Int row, column;
    if (description.field == kMatrixMarketFieldComplex) {
      double real_value, imag_value;
      if (!ReadMatrixMarketCoordinateComplexEntry(
              description, file, &row, &column, &real_value, &imag_value)) {
        result.reset();
        return result;
      }
      if ((mask == kEntryMaskLowerTriangle && row < column) ||
          (mask == kEntryMaskUpperTriangle && row > column)) {
        continue;
      }
      if (skip_explicit_zeros && real_value == 0. && imag_value == 0.) {
        ++num_skipped_entries;
        continue;
      }
    } else if (description.field == kMatrixMarketFieldReal) {
      double value;
      if (!ReadMatrixMarketCoordinateRealEntry(description, file, &row, &column,
                                               &value)) {
        result.reset();
        return result;
      }
      if ((mask == kEntryMaskLowerTriangle && row < column) ||
          (mask == kEntryMaskUpperTriangle && row > column)) {
        continue;
      }
      if (skip_explicit_zeros && value == 0.) {
        ++num_skipped_entries;
        continue;
      }
    } else {
      if (!ReadMatrixMarketCoordinateIndices(description, file, &row,
                                             &column)) {
        result.reset();
        return result;
      }
      if ((mask == kEntryMaskLowerTriangle && row < column) ||
          (mask == kEntryMaskUpperTriangle && row > column)) {
        continue;
      }
    }

    result->QueueEdgeAddition(row, column);
    if (row != column && description.symmetry != kMatrixMarketSymmetryGeneral) {
      result->QueueEdgeAddition(column, row);
    }
  }
  result->FlushEdgeQueues();

  if (num_skipped_entries) {
    std::cout << "Skipped " << num_skipped_entries
              << " explicitly zero entries." << std::endl;
  }

  return result;
}

inline void CoordinateGraph::ToMatrixMarket(const std::string& filename) const {
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Could not open " << filename << std::endl;
    return;
  }

  // Write the header.
  {
    std::ostringstream os;
    os << kMatrixMarketStampString << " " << kMatrixMarketObjectMatrixString
       << " " << kMatrixMarketFormatCoordinateString << " "
       << kMatrixMarketFieldRealString << " "
       << kMatrixMarketSymmetryGeneralString << "\n";
    file << os.str();
  }

  // Write the size information.
  {
    std::ostringstream os;
    os << NumSources() << " " << NumTargets() << " " << NumEdges() << "\n";
    file << os.str();
  }

  // Write out the entries.
  const Buffer<GraphEdge>& edges = Edges();
  for (const GraphEdge& edge : edges) {
    std::ostringstream os;
    // We must convert from 0-based to 1-based indexing.
    os << edge.first + 1 << " " << edge.second + 1 << " 1.\n";
    file << os.str();
  }
}

inline CoordinateGraph::~CoordinateGraph() {}

inline Int CoordinateGraph::NumSources() const QUOTIENT_NOEXCEPT {
  return num_sources_;
}

inline Int CoordinateGraph::NumTargets() const QUOTIENT_NOEXCEPT {
  return num_targets_;
}

inline Int CoordinateGraph::NumEdges() const QUOTIENT_NOEXCEPT {
  return edges_.Size();
}

inline void CoordinateGraph::Empty() {
  edges_.Clear();
  source_edge_offsets_.Clear();

  SwapClearVector(&edges_to_add_);
  SwapClearVector(&edges_to_remove_);

  num_sources_ = num_targets_ = 0;

  // Create a trivial source offset vector.
  source_edge_offsets_.Resize(1);
  source_edge_offsets_[0] = 0;
}

inline void CoordinateGraph::Resize(Int num_vertices) {
  AsymmetricResize(num_vertices, num_vertices);
}

inline void CoordinateGraph::AsymmetricResize(Int num_sources,
                                              Int num_targets) {
  if (num_sources == num_sources_ && num_targets == num_targets_) {
    return;
  }

  num_sources_ = num_sources;
  num_targets_ = num_targets;

  edges_.Clear();
  SwapClearVector(&edges_to_add_);
  SwapClearVector(&edges_to_remove_);

  source_edge_offsets_.Resize(num_sources + 1);
  for (Int source = 0; source <= num_sources; ++source) {
    source_edge_offsets_[source] = 0;
  }
}

inline void CoordinateGraph::ReserveEdgeAdditions(Int max_edge_additions) {
  edges_to_add_.reserve(max_edge_additions);
}

inline void CoordinateGraph::QueueEdgeAddition(Int source, Int target) {
  QUOTIENT_ASSERT(edges_to_add_.size() != edges_to_add_.capacity(),
                  "WARNING: Pushing back without first reserving space.");
  QUOTIENT_ASSERT(source >= 0 && source < num_sources_,
                  "ERROR: Source index was out of bounds.");
  QUOTIENT_ASSERT(target >= 0 && target < num_targets_,
                  "ERROR: Target index was out of bounds.");
  edges_to_add_.emplace_back(source, target);
}

inline void CoordinateGraph::FlushEdgeAdditionQueue(
    bool update_source_edge_offsets) {
  if (!edges_to_add_.empty()) {
    // Sort and erase duplicates from the list of edges to be added.
    std::sort(edges_to_add_.begin(), edges_to_add_.end());
    EraseDuplicatesInSortedVector(&edges_to_add_);

    // Perform a merge sort.
    const Buffer<GraphEdge> edges_copy(edges_);
    edges_.Resize(edges_copy.Size() + edges_to_add_.size());
    std::merge(edges_copy.begin(), edges_copy.end(), edges_to_add_.begin(),
               edges_to_add_.end(), edges_.begin());
    SwapClearVector(&edges_to_add_);

    // Erase duplicate edges.
    EraseDuplicatesInSortedVector(&edges_);
  }

  if (update_source_edge_offsets) {
    UpdateSourceEdgeOffsets();
  }
}

inline void CoordinateGraph::ReserveEdgeRemovals(Int max_edge_removals) {
  edges_to_remove_.reserve(max_edge_removals);
}

inline void CoordinateGraph::QueueEdgeRemoval(Int source, Int target) {
  QUOTIENT_ASSERT(source >= 0 && source < num_sources_,
                  "ERROR: Source index was out of bounds.");
  QUOTIENT_ASSERT(target >= 0 && target < num_targets_,
                  "ERROR: Target index was out of bounds.");
  edges_to_remove_.emplace_back(source, target);
}

inline void CoordinateGraph::FlushEdgeRemovalQueue(
    bool update_source_edge_offsets) {
  if (!edges_to_remove_.empty()) {
    // Sort and erase duplicates from the list of edges to be removed.
    std::sort(edges_to_remove_.begin(), edges_to_remove_.end());
    EraseDuplicatesInSortedVector(&edges_to_remove_);

    const Int num_edges = edges_.Size();
    Int num_packed = 0;
    for (Int index = 0; index < num_edges; ++index) {
      auto iter = std::lower_bound(edges_to_remove_.begin(),
                                   edges_to_remove_.end(), edges_[index]);
      if (iter == edges_to_remove_.end() || *iter != edges_[index]) {
        // The current edge should be kept, so pack it from the left.
        edges_[num_packed++] = edges_[index];
      }
    }
    edges_.Resize(num_packed);
    SwapClearVector(&edges_to_remove_);
  }

  if (update_source_edge_offsets) {
    UpdateSourceEdgeOffsets();
  }
}

inline void CoordinateGraph::FlushEdgeQueues() {
  if (EdgeQueuesAreEmpty()) {
    // Skip the recomputation of the source offsets.
    return;
  }
  FlushEdgeRemovalQueue(false /* update_source_edge_offsets */);
  FlushEdgeAdditionQueue(true /* update_source_edge_offsets */);
}

inline bool CoordinateGraph::EdgeQueuesAreEmpty() const QUOTIENT_NOEXCEPT {
  return edges_to_add_.empty() && edges_to_remove_.empty();
}

inline void CoordinateGraph::AddEdge(Int source, Int target) {
  ReserveEdgeAdditions(1);
  QueueEdgeAddition(source, target);
  FlushEdgeQueues();
}

inline void CoordinateGraph::RemoveEdge(Int source, Int target) {
  ReserveEdgeRemovals(1);
  QueueEdgeRemoval(source, target);
  FlushEdgeQueues();
}

inline const GraphEdge& CoordinateGraph::Edge(Int edge_index) const
    QUOTIENT_NOEXCEPT {
  return edges_[edge_index];
}

inline const Buffer<GraphEdge>& CoordinateGraph::Edges() const
    QUOTIENT_NOEXCEPT {
  return edges_;
}

inline Int CoordinateGraph::SourceEdgeOffset(Int source) const
    QUOTIENT_NOEXCEPT {
  QUOTIENT_ASSERT(
      EdgeQueuesAreEmpty(),
      "Tried to retrieve a source edge offset when edge queues weren't empty");
  return source_edge_offsets_[source];
}

inline Int CoordinateGraph::EdgeOffset(Int source,
                                       Int target) const QUOTIENT_NOEXCEPT {
  const Int source_edge_offset = SourceEdgeOffset(source);
  const Int next_source_edge_offset = SourceEdgeOffset(source + 1);
  auto iter = std::lower_bound(edges_.begin() + source_edge_offset,
                               edges_.begin() + next_source_edge_offset,
                               GraphEdge(source, target));
  return std::distance(edges_.begin(), iter);
}

inline bool CoordinateGraph::EdgeExists(Int source,
                                        Int target) const QUOTIENT_NOEXCEPT {
  const Int index = EdgeOffset(source, target);
  const GraphEdge& edge = Edge(index);
  return edge.first == source && edge.second == target;
}

inline Int CoordinateGraph::NumConnections(Int source) const QUOTIENT_NOEXCEPT {
  return SourceEdgeOffset(source + 1) - SourceEdgeOffset(source);
}

inline void CoordinateGraph::UpdateSourceEdgeOffsets() {
  const Int num_edges = edges_.Size();
  source_edge_offsets_.Resize(num_sources_ + 1);
  Int source_edge_offset = 0;
  Int prev_source = -1;
  for (Int edge_index = 0; edge_index < num_edges; ++edge_index) {
    const Int source = edges_[edge_index].first;
    QUOTIENT_ASSERT(source >= prev_source, "Sources were not properly sorted.");

    // Fill in the source offsets from prev_source to source - 1.
    for (; prev_source < source; ++prev_source) {
      source_edge_offsets_[source_edge_offset++] = edge_index;
    }
  }

  // Fill in the end of the source offset buffer.
  for (; source_edge_offset <= num_sources_; ++source_edge_offset) {
    source_edge_offsets_[source_edge_offset] = num_edges;
  }
}

}  // namespace quotient

#endif  // ifndef QUOTIENT_COORDINATE_GRAPH_IMPL_H_
