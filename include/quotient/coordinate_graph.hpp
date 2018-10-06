/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_COORDINATE_GRAPH_H_
#define QUOTIENT_COORDINATE_GRAPH_H_

#include <algorithm>
#include <fstream>
#include <memory>
#include <iostream>
#include <sstream>
#include <vector>

#include "quotient/config.hpp"

namespace quotient {

// A pairing of the source and target vertices of a graph edge.
typedef std::pair<Int, Int> GraphEdge;

// The character that Matrix Market comment lines begin with.
static constexpr char kMatrixMarketCommentChar = '%';

// The Matrix Market 'complex' field string.
static constexpr char kMatrixMarketFieldComplexString[] = "complex";

// The Matrix Market 'double' field string.
static constexpr char kMatrixMarketFieldDoubleString[] = "real";

// The Matrix Market 'integer' field string.
static constexpr char kMatrixMarketFieldIntegerString[] = "integer";

// The Matrix Market 'pattern' field string.
static constexpr char kMatrixMarketFieldPatternString[] = "pattern";

// The Matrix Market 'real' field string.
static constexpr char kMatrixMarketFieldRealString[] = "real";

// The Matrix Market 'array' format string.
static constexpr char kMatrixMarketFormatArrayString[]= "array";

// The Matrix Market 'coordinate' format string.
static constexpr char kMatrixMarketFormatCoordinateString[] = "coordinate";

// The Matrix Market 'matrix' object string.
static constexpr char kMatrixMarketObjectMatrixString[] = "matrix";

// The Matrix Market 'vector' object string.
static constexpr char kMatrixMarketObjectVectorString[] = "vector";

// The Matrix Market 'hermitian' symmetry string.
static constexpr char kMatrixMarketSymmetryHermitianString[] = "hermitian";

// The Matrix Market 'general' symmetry string.
static constexpr char kMatrixMarketSymmetryGeneralString[] = "general";

// The Matrix Market 'skew-symmetric' symmetry string.
static constexpr char kMatrixMarketSymmetrySkewSymmetricString[] =
    "skew-symmetric";

// The Matrix Market 'symmetric' symmetry string.
static constexpr char kMatrixMarketSymmetrySymmetricString[] = "symmetric";

// The Matrix Market stamp string.
static constexpr char kMatrixMarketStampString[] = "%%MatrixMarket";


// A helper routine for freeing the memory of an std::vector.
template<typename T>
void SwapClearVector(std::vector<T>* vec) {
  std::vector<T>().swap(*vec);
}


// Removes all duplicate entries from a sorted vector.
template<typename T>
void EraseDuplicatesInSortedVector(std::vector<T>* vec) {
  vec->erase(std::unique(vec->begin(), vec->end()), vec->end());
}


// A coordinate-format graph data structure that supports different source
// and target set sizes. The primary storage is a lexicographically sorted
// std::vector<GraphEdge> and an associated std::vector<Int> of source offsets
// (which serve the same role as in Compressed Sparse Row (CSR) format). Thus,
// this storage scheme is a superset of the CSR format that explicitly stores
// both source and target indices for each edge.
//
// The class is designed so that the sorting and offset computation overhead
// can be amortized over batches of edge additions and removals.
//
// For example, the code block:
//
//   quotient::CoordinateGraph graph;
//   graph.Resize(5);
//   graph.ReserveEdgeAdditions(6);
//   graph.QueueEdgeAddition(3, 4);
//   graph.QueueEdgeAddition(2, 3);
//   graph.QueueEdgeAddition(2, 0);
//   graph.QueueEdgeAddition(4, 2);
//   graph.QueueEdgeAddition(4, 4);
//   graph.QueueEdgeAddition(3, 2);
//   graph.FlushEdgeQueues();
//   const std::vector<quotient::GraphEdge>& edges = graph.Edges();
//
// would return a reference to the underlying std::vector<quotient::GraphEdge>
// of 'graph', which should contain the edge sequence:
//   (2, 0), (2, 3), (3, 2), (3, 4), (4, 2), (4, 4).
//
// Similarly, subsequently running the code block:
//
//   graph.ReserveEdgeRemovals(2);
//   graph.QueueEdgeRemoval(2, 3);
//   graph.QueueEdgeRemoval(0, 4);
//   graph.FlushEdgeQueues();
//
// would modify the std::vector underlying the 'edges' reference to now
// contain the edge sequence:
//   (2, 0), (3, 2), (3, 4), (4, 2), (4, 4).
//
// TODO(Jack Poulson): Add support for 'END' index marker so that ranges
// can be easily incorporated.
class CoordinateGraph {
 public:
  // The trivial constructor.
  CoordinateGraph();

  // The copy constructor.
  CoordinateGraph(const CoordinateGraph& graph);

  // The assignment operator.
  const CoordinateGraph& operator=(const CoordinateGraph& graph);

  // Builds and returns a CoordinateGraph from a Matrix Market description.
  static std::unique_ptr<CoordinateGraph> FromMatrixMarket(
      const std::string& filename);

  // A trivial destructor.
  ~CoordinateGraph();

  // Returns the ground set size of the source vertices.
  Int NumSources() const noexcept;

  // Returns the ground set size of the target vertices.
  Int NumTargets() const noexcept;

  // Returns the number of edges in the graph.
  Int NumEdges() const noexcept;

  // Removes all edges and changes the source and target vertex ground set
  // sizes to zero.
  void Empty(bool free_resources);

  // Changes both the source and target ground set sizes to 'num_vertices'.
  void Resize(Int num_vertices);

  // Chnages the source ground set size to 'num_sources', and likewise for the
  // target ground set size and 'num_targets'.
  void AsymmetricResize(Int num_sources, Int num_targets);

  // Allocates space so that up to 'max_edge_additions' calls to
  // 'QueueEdgeAddition' can be performed without another memory allocation.
  void ReserveEdgeAdditions(Int max_edge_additions);

  // Appends the edge (source, target) to the edge list without putting the
  // edge list in lexicographic order or updating the source offsets.
  void QueueEdgeAddition(Int source, Int target);

  // Allocates space so that up to 'max_edge_additions' calls to
  // 'QueueEdgeRemoval' can be performed without another memory allocation.
  void ReserveEdgeRemovals(Int max_edge_removals);

  // Appends the edge (source, target) to the list of edges to be removed.
  void QueueEdgeRemoval(Int source, Int target);

  // All queued edge additions and removals are applied, the edge list is
  // lexicographically sorted, and the source offsets are then updated.
  void FlushEdgeQueues();

  // Returns true if there are no edges queued for addition or removal.
  bool EdgeQueuesAreEmpty() const noexcept;

  // Adds the edge (source, target) into the graph.
  //
  // NOTE: This routine involves a merge sort involving all of the edges. It is
  // preferable to amortize this cost by batching together several edge
  // additions.
  void AddEdge(Int source, Int target);

  // Removes the edge (source, target) into the graph.
  //
  // NOTE: This routine can involve a merge sort involving all of the edges.
  // It is preferable to amortize this cost by batching together several edge
  // removals.
  void RemoveEdge(Int source, Int target);

  // Returns a reference to the edge with the given index.
  const GraphEdge& Edge(Int edge_index) const noexcept;

  // Returns a reference to the underlying vector of edges.
  const std::vector<GraphEdge>& Edges() const noexcept;

  // Returns the offset into the edge vector where edges from the given source
  // begin.
  Int SourceEdgeOffset(Int source) const noexcept;

  // Returns the offset into the edge vector where the (source, target) edge
  // would be inserted.
  Int EdgeOffset(Int source, Int target) const noexcept;

  // Returns true if the (source, target) edge exists.
  bool EdgeExists(Int source, Int target) const noexcept;

  // Returns the number of targets that the given source is connected to.
  Int NumConnections(Int source) const noexcept;

 private:

  enum MatrixMarketObject {
    kMatrixMarketObjectMatrix,
    kMatrixMarketObjectVector,
  };

  enum MatrixMarketFormat {
    kMatrixMarketFormatArray,
    kMatrixMarketFormatCoordinate,
  };

  enum MatrixMarketField {
    kMatrixMarketFieldReal,
    kMatrixMarketFieldComplex,
    kMatrixMarketFieldInteger,
    kMatrixMarketFieldPattern,
  };

  enum MatrixMarketSymmetry {
    kMatrixMarketSymmetryGeneral,
    kMatrixMarketSymmetrySymmetric,
    kMatrixMarketSymmetrySkewSymmetric,
    kMatrixMarketSymmetryHermitian,
  };

  struct MatrixMarketDescription {
    MatrixMarketObject object;

    MatrixMarketFormat format;

    MatrixMarketField field;

    MatrixMarketSymmetry symmetry;

    // A trivial constructor.
    MatrixMarketDescription();

    // Builds the MatrixMarketDescription by parsing the header line of the
    // Matrix Market file. Returns true if the parse was successful.
    bool ParseFromHeaderLine(const std::string& header_line);
  };

  // The ground set size for the source vertices.
  Int num_sources_;

  // The ground set size for the target vertices.
  Int num_targets_;

  // The (lexicographically sorted) list of edges in the graph.
  std::vector<GraphEdge> edges_;

  // A list of length num_sources_ + 1, where source_edge_offsets_[source]
  // indicates the location in edges_ where the edge (source, 0) would be
  // inserted.
  std::vector<Int> source_edge_offsets_;

  // The list of edges currently queued for addition into the graph.
  std::vector<GraphEdge> edges_to_add_;

  // The list of edges currently queued for removal from the graph.
  std::vector<GraphEdge> edges_to_remove_;

  // Incorporates the edges currently residing in 'edges_to_add_' into the
  // graph (and then clears 'edges_to_add_').
  //
  // If 'update_source_edge_offsets' is true, then 'source_edge_offsets_' is
  // recomputed.
  void FlushEdgeAdditionQueue(bool update_source_edge_offsets);

  // Removes the edges residing in 'edges_to_remove_' from the graph (and then
  // clears 'edges_to_remove_').
  //
  // If 'update_source_edge_offsets' is true, then 'source_edge_offsets_' is
  // recomputed.
  void FlushEdgeRemovalQueue(bool update_source_edge_offsets);

  // Recomputes 'source_edge_offsets_' based upon the current value of 'edges_'.
  void UpdateSourceEdgeOffsets();
};

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
  edges_to_add_ = graph.edges_to_add_;
  edges_to_remove_ = graph.edges_to_remove_;

  return *this;
}

CoordinateGraph::MatrixMarketDescription::MatrixMarketDescription() { }

bool CoordinateGraph::MatrixMarketDescription::ParseFromHeaderLine(
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

std::unique_ptr<CoordinateGraph> CoordinateGraph::FromMatrixMarket(
      const std::string& filename) {
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
    result->QueueEdgeAddition(source, target);
    if (source != target &&
        description.symmetry != kMatrixMarketSymmetryGeneral) {
      result->QueueEdgeAddition(target, source);
    }
  }
  result->FlushEdgeQueues();

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

#endif // ifndef QUOTIENT_COORDINATE_GRAPH_H_
