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
#include <memory>
#include <iostream>
#include <sstream>
#include <vector>

#include "quotient/config.hpp"
#include "quotient/coordinate_graph.hpp"

namespace quotient {

template<typename T>
void SwapClearVector(std::vector<T>* vec) {
  std::vector<T>().swap(*vec);
}

template<typename T>
void EraseDuplicatesInSortedVector(std::vector<T>* vec) {
  vec->erase(std::unique(vec->begin(), vec->end()), vec->end());
}

inline CoordinateGraph::CoordinateGraph() 
: num_sources_(0), num_targets_(0) { }

inline CoordinateGraph::CoordinateGraph(const CoordinateGraph& graph) {
  if (&graph == this) {
    return;
  }
  *this = graph;
}

inline const CoordinateGraph& CoordinateGraph::operator=(
    const CoordinateGraph& graph) {
  if (&graph == this) {
    return *this;
  }

  num_sources_ = graph.num_sources_;
  num_targets_ = graph.num_targets_;
  edges_ = graph.edges_;
  source_edge_offsets_ = graph.source_edge_offsets_;
  edges_to_add_ = graph.edges_to_add_;
  edges_to_remove_ = graph.edges_to_remove_;

  return *this;
}

inline CoordinateGraph::MatrixMarketDescription::MatrixMarketDescription() { }

inline bool CoordinateGraph::MatrixMarketDescription::ParseFromHeaderLine(
    const std::string& header_line) {
  std::stringstream line_stream(header_line);

  // Extract the tokens of the header line.
  std::string stamp;
  line_stream >> stamp;
  if (stamp != kMatrixMarketStampString) {
    std::cerr << "Invalid Matrix Market stamp." << std::endl;
    return false;
  }
  std::string object_string;
  if (!(line_stream >> object_string)) {
    std::cerr << "Missing Matrix Market object." << std::endl;
    return false;
  }
  std::string format_string;
  if (!(line_stream >> format_string)) {
    std::cerr << "Missing Matrix Market format." << std::endl;
    return false;
  }
  std::string field_string;
  if (!(line_stream >> field_string)) {
    std::cerr << "Missing Matrix Market field." << std::endl;
    return false;
  }
  std::string symmetry_string;
  if (!(line_stream >> symmetry_string)) {
    std::cerr << "Missing Matrix Market symmetry." << std::endl;
    return false;
  }

  // Determine the type of object.
  if (object_string == kMatrixMarketObjectMatrixString) { 
    object = kMatrixMarketObjectMatrix;
  } else if (object_string == kMatrixMarketObjectVectorString) {
    object = kMatrixMarketObjectVector;
  } else {
    std::cerr << "Invalid Matrix Market object string." << std::endl;
    return false;
  }

  // Determine the storage format.
  if (format_string == kMatrixMarketFormatArrayString) {
    format = kMatrixMarketFormatArray;
  } else if (format_string == kMatrixMarketFormatCoordinateString) {
    format = kMatrixMarketFormatCoordinate;
  } else {
    std::cerr << "Invalid Matrix Market format string." << std::endl;
    return false;
  }

  // Determine the underlying field.
  if (field_string == kMatrixMarketFieldComplexString) {
    field = kMatrixMarketFieldComplex;
  } else if (field_string == kMatrixMarketFieldDoubleString ||
             field_string == kMatrixMarketFieldRealString) {
    field = kMatrixMarketFieldReal;
  } else if (field_string == kMatrixMarketFieldIntegerString) {
    field = kMatrixMarketFieldInteger;
  } else if (field_string == kMatrixMarketFieldPatternString) {
    field = kMatrixMarketFieldPattern;
  } else {
    std::cerr << "Invalid Matrix Market field string." << std::endl;
    return false;
  }

  // Determine the symmetry.
  if (symmetry_string == kMatrixMarketSymmetryGeneralString) {
    symmetry = kMatrixMarketSymmetryGeneral;
  } else if (symmetry_string == kMatrixMarketSymmetrySymmetricString) {
    symmetry = kMatrixMarketSymmetrySymmetric;
  } else if (symmetry_string == kMatrixMarketSymmetrySkewSymmetricString) {
    symmetry = kMatrixMarketSymmetrySkewSymmetric;
  } else if (symmetry_string == kMatrixMarketSymmetryHermitianString) {
    symmetry = kMatrixMarketSymmetryHermitian;
  } else {
    std::cerr << "Invalid Matrix Market symmetry string." << std::endl;
    return false;
  }

  // Validate the description.
  if (format == kMatrixMarketFormatArray &&
      field == kMatrixMarketFieldPattern) {
    std::cerr << "The 'array' format and 'pattern' field are incompatible."
              << std::endl;
    return false;
  }
  if (field == kMatrixMarketFieldPattern &&
      symmetry == kMatrixMarketSymmetrySkewSymmetric) {
    // Due to the line
    //
    //   "symmetry is either general (legal for real, complex, integer or
    //    pattern fields), symmetric (real, complex, integer or pattern),
    //    skew-symmetric (real, complex or integer), or hermitian (complex
    //    only)."
    //
    // from http://people.sc.fsu.edu/~jburkardt/data/mm/mm.html, we disallow
    // patterned skew-symmetry.
    std::cerr << "The 'pattern' field and 'skew-symmetric' symmetry are "
                 "incompatible." << std::endl;
    return false;
  }
  if (field != kMatrixMarketFieldComplex &&
      symmetry == kMatrixMarketSymmetryHermitian) {
    std::cerr << "The hermitian symmetry technically requires complex data."
              << std::endl;
  }

  return true;
}

inline std::unique_ptr<CoordinateGraph> CoordinateGraph::FromMatrixMarket(
    const std::string& filename,
    bool skip_explicit_zeros,
    EntryMask mask) {
  std::unique_ptr<CoordinateGraph> result;
  std::ifstream file(filename);
  if (!file.is_open()) { 
    std::cerr << "Could not open " << filename << std::endl;
    return result;
  }
  
  // Fill the description of the Matrix Market data.
  std::string line;
  if (!std::getline(file, line)) {
    std::cerr << "Could not read header line from Matrix Market file."
              << std::endl;
    return result;
  }
  MatrixMarketDescription description;
  if (!description.ParseFromHeaderLine(line)) {
    std::cerr << "Could not parse header line from Matrix Market file."
              << std::endl;
    return result;
  }
  if (description.object == kMatrixMarketObjectVector) {
    std::cerr << "The Matrix Market 'vector' object is incompatible with "
                 "CoordinateGraph." << std::endl;
    return result;
  }

  // Skip the comment lines.
  while (file.peek() == kMatrixMarketCommentChar) {
    std::getline(file, line);
  }
 
  result.reset(new CoordinateGraph);
  if (description.format == kMatrixMarketFormatArray) {
    // Read the size of the matrix.
    Int matrix_height, matrix_width;
    if (!std::getline(file, line)) {
      std::cerr << "Could not extract the array size line." << std::endl;
      result.reset();
      return result;
    }
    std::stringstream line_stream(line);
    if (!(line_stream >> matrix_height)) {  
      std::cerr << "Missing matrix height in Matrix Market file." << std::endl;
      result.reset();
      return result;
    }
    if (description.object == kMatrixMarketObjectMatrix) {
      if (!(line_stream >> matrix_width)) {
        std::cerr << "Missing matrix width in Matrix Market file." << std::endl;
        result.reset();
        return result;
      }
    } else {
      matrix_width = 1;
    }
    
    // Fill a fully-connected graph.
    result->AsymmetricResize(matrix_height, matrix_width);
    result->ReserveEdgeAdditions(matrix_height * matrix_width);
    for (Int source = 0; source < matrix_height; ++source) {
      for (Int target = 0; target < matrix_width; ++target) {
        result->QueueEdgeAddition(source, target);
      }
    }
    result->FlushEdgeQueues();
    return result;
  }

  // Read in the number of matrix dimensions and the number of explicit
  // nonzeros specified in the file.
  Int matrix_height, matrix_width, num_explicit_nonzeros;
  if (!std::getline(file, line)) {
    std::cerr << "Could not extract the coordinate size line." << std::endl;
    result.reset();
    return result;
  }
  std::stringstream line_stream(line); 
  if (!(line_stream >> matrix_height)) {
    std::cerr << "Missing matrix height in Matrix Market file." << std::endl;
    result.reset();
    return result;
  }
  if (description.object == kMatrixMarketObjectMatrix) {
    if (!(line_stream >> matrix_width)) {
      std::cerr << "Missing matrix width in Matrix Market file." << std::endl;
      result.reset();
      return result;
    }
  } else {
    matrix_width = 1;
  }
  if (!(line_stream >> num_explicit_nonzeros)) {
    std::cerr << "Missing num_nonzeros in Matrix Market file." << std::endl;
    result.reset();
    return result;
  }

  // Fill in the edges.
  Int num_skipped_edges = 0;
  const Int num_nonzeros_bound =
      description.symmetry == kMatrixMarketSymmetryGeneral ?
      num_explicit_nonzeros : 2 * num_explicit_nonzeros;
  result->AsymmetricResize(matrix_height, matrix_width);
  result->ReserveEdgeAdditions(num_nonzeros_bound);
  for (Int edge_index = 0; edge_index < num_explicit_nonzeros; ++edge_index) {
    Int source, target;
    if (!std::getline(file, line)) {
      std::cerr << "Could not extract nonzero description from Matrix Market "
                   "file." << std::endl;
      result.reset();
      return result;
    }
    std::stringstream line_stream(line);
    if (!(line_stream >> source)) {
      std::cerr << "Could not extract row index of nonzero." << std::endl;
      result.reset();
      return result;
    }
    --source; // Convert from 1-based to 0-based indexing.
    if (description.object == kMatrixMarketObjectMatrix) {
      if (!(line_stream >> target)) {
        std::cerr << "Could not extract column index of nonzero." << std::endl;
        result.reset();
        return result;
      }
      --target; // Convert from 1-based to 0-based indexing.
    } else {
      target = 0;
    }
    if ((mask == kEntryMaskLowerTriangle && source < target) ||
        (mask == kEntryMaskUpperTriangle && source > target)) {
      continue;
    }

    if (skip_explicit_zeros) {
      // Skip this entry if it is numerically zero.
      if (description.field == kMatrixMarketFieldReal) {
        double value;
        if (!(line_stream >> value)) {
          std::cerr << "Could not extract real value of nonzero." << std::endl;
        }
        if (value == 0.) {
          ++num_skipped_edges;
          continue;
        }
      } else if (description.field == kMatrixMarketFieldComplex) {
        double real_value, imag_value;
        if (!(line_stream >> real_value)) {
          std::cerr << "Could not extract real value of nonzero." << std::endl;
        }
        if (!(line_stream >> imag_value)) {
          std::cerr << "Could not extract imag value of nonzero." << std::endl;
        }
        if (real_value == 0. && imag_value == 0.) {
          ++num_skipped_edges;
          continue;
        }
      }
    }

    result->QueueEdgeAddition(source, target);
    if (source != target &&
        description.symmetry != kMatrixMarketSymmetryGeneral) {
      result->QueueEdgeAddition(target, source);
    }
  }
  result->FlushEdgeQueues();

  if (skip_explicit_zeros && num_skipped_edges > 0) {
    std::cout << "Skipped " << num_skipped_edges << " explicitly zero edges."
              << std::endl;
  }

  return result;
}

inline CoordinateGraph::~CoordinateGraph() { }

inline Int CoordinateGraph::NumSources() const noexcept { return num_sources_; }

inline Int CoordinateGraph::NumTargets() const noexcept { return num_targets_; }

inline Int CoordinateGraph::NumEdges() const noexcept { return edges_.size(); }

inline void CoordinateGraph::Empty(bool free_resources) {
  if (free_resources) {
    SwapClearVector(&edges_);
    SwapClearVector(&source_edge_offsets_);
    SwapClearVector(&edges_to_add_);
    SwapClearVector(&edges_to_remove_);
  } else {
    edges_.clear();
    edges_to_add_.clear();
    edges_to_remove_.clear();
  }

  num_sources_ = num_targets_ = 0;

  // Create a trivial source offset vector.
  source_edge_offsets_.resize(1);
  source_edge_offsets_[0] = 0;
}

inline void CoordinateGraph::Resize(Int num_vertices) {
  AsymmetricResize(num_vertices, num_vertices);
}

inline void CoordinateGraph::AsymmetricResize(
    Int num_sources, Int num_targets) {
  if (num_sources == num_sources_ && num_targets == num_targets_) {
    return;
  }

  num_sources_ = num_sources;
  num_targets_ = num_targets;

  edges_.clear();
  edges_to_add_.clear();
  edges_to_remove_.clear();

  source_edge_offsets_.resize(num_sources + 1);
  for (Int source = 0; source <= num_sources; ++source) {
    source_edge_offsets_[source] = 0;
  }
}

inline void CoordinateGraph::ReserveEdgeAdditions(Int max_edge_additions) {
  edges_to_add_.reserve(max_edge_additions);
}

inline void CoordinateGraph::QueueEdgeAddition(Int source, Int target) {
#ifdef QUOTIENT_DEBUG
  if (edges_to_add_.size() == edges_to_add_.capacity()) {  
    std::cerr << "WARNING: Pushing back without first reserving space."
	      << std::endl;
  }
  if (source < 0 || source >= num_sources_) {
    std::cerr << "ERROR: Source index was out of bounds." << std::endl;
    return;
  }
  if (target < 0 || target >= num_targets_) {
    std::cerr << "ERROR: Target index was out of bounds." << std::endl;
    return;
  }
#endif // ifdef QUOTIENT_DEBUG
  edges_to_add_.emplace_back(source, target);
}

inline void CoordinateGraph::FlushEdgeAdditionQueue(
    bool update_source_edge_offsets) {
  if (!edges_to_add_.empty()) {
    // Sort and erase duplicates from the list of edges to be added.
    std::sort(edges_to_add_.begin(), edges_to_add_.end());
    EraseDuplicatesInSortedVector(&edges_to_add_);

    // Perform a merge sort.
    const std::vector<GraphEdge> edges_copy(edges_);
    edges_.resize(0);
    edges_.resize(edges_copy.size() + edges_to_add_.size());
    std::merge(
      edges_copy.begin(), edges_copy.end(),
      edges_to_add_.begin(), edges_to_add_.end(),
      edges_.begin());
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
#ifdef QUOTIENT_DEBUG
  if (source < 0 || source >= num_sources_) {
    std::cerr << "ERROR: Source index was out of bounds." << std::endl;
    return;
  }
  if (target < 0 || target >= num_targets_) {
    std::cerr << "ERROR: Target index was out of bounds." << std::endl;
    return;
  }
#endif
  edges_to_remove_.emplace_back(source, target);
}

inline void CoordinateGraph::FlushEdgeRemovalQueue(
    bool update_source_edge_offsets) {
  if (!edges_to_remove_.empty()) {
    // Sort and erase duplicates from the list of edges to be removed.
    std::sort(edges_to_remove_.begin(), edges_to_remove_.end());
    EraseDuplicatesInSortedVector(&edges_to_remove_);

    const Int num_edges = edges_.size();
    Int num_removed = 0;
    for (Int index = 0; index < num_edges; ++index) { 
      auto iter = std::lower_bound(
        edges_to_remove_.begin(), edges_to_remove_.end(), edges_[index]);
      if (iter == edges_to_remove_.end() || *iter != edges_[index]) {
        // The current edge should be kept, so pack it from the left.
        edges_[index - num_removed] = edges_[index];
      } else {
        // The current edge should be skipped.
        ++num_removed;
      }
    }
    edges_.resize(num_edges - num_removed);
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

inline bool CoordinateGraph::EdgeQueuesAreEmpty() const noexcept {
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

inline const GraphEdge& CoordinateGraph::Edge(Int edge_index) const noexcept {
#ifdef QUOTIENT_DEBUG
  return edges_.at(edge_index);
#else
  return edges_[edge_index];
#endif
}

inline const std::vector<GraphEdge>& CoordinateGraph::Edges() const noexcept {
  return edges_;
}

inline Int CoordinateGraph::SourceEdgeOffset(Int source) const noexcept {
#ifdef QUOTIENT_DEBUG
  if (!EdgeQueuesAreEmpty()) {
    std::cerr << "Tried to retrieve a source edge offset when edge queues "
                 "were not empty." << std::endl;
    return -1;
  }
  return source_edge_offsets_.at(source);
#else
  return source_edge_offsets_[source];
#endif
}

inline Int CoordinateGraph::EdgeOffset(Int source, Int target) const noexcept {
  const Int source_edge_offset = SourceEdgeOffset(source);
  const Int next_source_edge_offset = SourceEdgeOffset(source + 1);
  auto iter = std::lower_bound(
      edges_.begin() + source_edge_offset,
      edges_.begin() + next_source_edge_offset,
      GraphEdge(source, target));
  return iter - edges_.begin();
}

inline bool CoordinateGraph::EdgeExists(Int source, Int target) const noexcept {
  const Int index = EdgeOffset(source, target);
  const GraphEdge& edge = Edge(index);
  return edge.first == source && edge.second == target;
}

inline Int CoordinateGraph::NumConnections(Int source) const noexcept {
  return SourceEdgeOffset(source + 1) - SourceEdgeOffset(source);
}

inline void CoordinateGraph::UpdateSourceEdgeOffsets() {
  const Int num_edges = edges_.size();
  source_edge_offsets_.resize(num_sources_ + 1);
  Int source_edge_offset = 0;  
  Int prev_source = -1;
  for (Int edge_index = 0; edge_index < num_edges; ++edge_index) {
    const Int source = edges_[edge_index].first;
#ifdef QUOTIENT_DEBUG
    if (source < prev_source) {
      std::cerr << "Sources were not properly sorted." << std::endl;
      return;
    }
#endif

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

} // namespace quotient

#endif // ifndef QUOTIENT_COORDINATE_GRAPH_IMPL_H_
