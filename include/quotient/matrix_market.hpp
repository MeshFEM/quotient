/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_MATRIX_MARKET_H_
#define QUOTIENT_MATRIX_MARKET_H_

#include "quotient/config.hpp"

namespace quotient {

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

} // namespace quotient

#include "quotient/matrix_market-impl.hpp"

#endif // ifndef QUOTIENT_MATRIX_MARKET_H_
