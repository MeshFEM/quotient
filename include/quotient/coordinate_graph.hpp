/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_COORDINATE_GRAPH_H_
#define QUOTIENT_COORDINATE_GRAPH_H_

#include <memory>
#include <vector>

#include "quotient/config.hpp"

namespace quotient {

// A pairing of the source and target vertices of a graph edge.
typedef std::pair<Int, Int> GraphEdge;

// An enum for representing a portion of a square matrix.
enum EntryMask {
  // Allow all entries.
  kEntryMaskFull,

  // Allow only the entries in the lower triangle.
  kEntryMaskLowerTriangle,

  // Allow only the entries in the upper triangle.
  kEntryMaskUpperTriangle,
};

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
void SwapClearVector(std::vector<T>* vec);


// Removes all duplicate entries from a sorted vector.
template<typename T>
void EraseDuplicatesInSortedVector(std::vector<T>* vec);


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
      const std::string& filename,
      bool skip_explicit_zeros,
      EntryMask mask=kEntryMaskFull);

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

  // A representation of the 'Object' options of a Matrix Market file.
  enum MatrixMarketObject {
    // The object represents a matrix and so both row and column indices should
    // be provided.
    kMatrixMarketObjectMatrix,

    // The object represents a vector and so only row indices should be
    // provided.
    kMatrixMarketObjectVector,
  };

  // A representation of the 'Format' options of a Matrix Market file.
  enum MatrixMarketFormat {
    // The matrix is treated as dense and individual indices will not be
    // provided.
    kMatrixMarketFormatArray,

    // The matrix is treated as sparse and individual indices will be
    // provided.
    kMatrixMarketFormatCoordinate,
  };

  // A representation of the 'Field' options of a Matrix Market file.
  enum MatrixMarketField {
    // The matrix contains real, double-precision entries.
    kMatrixMarketFieldReal,

    // The matrix contains complex, double-precision entries.
    kMatrixMarketFieldComplex,

    // The matrix contains integer entries.
    kMatrixMarketFieldInteger,

    // The matrix does not have explicitly specified numerical values.
    kMatrixMarketFieldPattern,
  };

  // A representation of the 'Symmetry' options of a Matrix Market file.
  enum MatrixMarketSymmetry {
    // No symmetry is assumed.
    kMatrixMarketSymmetryGeneral,

    // The matrix is assumed symmetric.
    kMatrixMarketSymmetrySymmetric,

    // The matrix is assumed skew-symmetric.
    kMatrixMarketSymmetrySkewSymmetric,

    // The matrix is assumed Hermitian.
    kMatrixMarketSymmetryHermitian,
  };

  // A representation of a Matrix Market file's metadata.
  struct MatrixMarketDescription {
    // Whether the file contains a matrix or a vector.
    MatrixMarketObject object;

    // Whether the corresponding matrix is assumed dense or sparse.
    MatrixMarketFormat format;

    // The type of numerical values associated with the matrix.
    MatrixMarketField field;

    // The assumed symmetry (if any) of the matrix.
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

} // namespace quotient

#include "quotient/coordinate_graph-impl.hpp"

#endif // ifndef QUOTIENT_COORDINATE_GRAPH_H_
