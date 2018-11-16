/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_MATRIX_MARKET_IMPL_H_
#define QUOTIENT_MATRIX_MARKET_IMPL_H_

#include <memory>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "quotient/matrix_market.hpp"

namespace quotient {

inline MatrixMarketDescription::MatrixMarketDescription() { }

inline bool MatrixMarketDescription::ParseFromHeaderLine(
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

inline bool ReadMatrixMarketDescription(
    std::ifstream& file, MatrixMarketDescription* description) {
  // Get the first line of the file.
  std::string line;
  if (!std::getline(file, line)) {
    std::cerr << "Could not read header line from Matrix Market file."
              << std::endl;
    return false;
  }

  // Parse the first line of the file into the description.
  if (!description->ParseFromHeaderLine(line)) {
    std::cerr << "Could not parse header line from Matrix Market file."
              << std::endl;
    return false;
  }

  // Perform a preliminary consistency check.
  if (description->object == kMatrixMarketObjectVector) {
    std::cerr << "The Matrix Market 'vector' object is incompatible with "
                 "CoordinateGraph." << std::endl;
    return false;
  }

  return true;
}

inline bool ReadMatrixMarketArrayMetadata(
    const MatrixMarketDescription& description,
    std::ifstream& file,
    Int* num_rows,
    Int* num_columns) {
  std::string line;

  // Skip any comment lines.
  while (file.peek() == kMatrixMarketCommentChar) {
    std::getline(file, line);
  }

  // Get a stringstream for the relevant line of the file.
  if (!std::getline(file, line)) {
    std::cerr << "Could not extract the array metadata line." << std::endl;
    return false;
  }
  std::stringstream line_stream(line);

  // Read the number of rows.
  if (!(line_stream >> *num_rows)) {
    std::cerr << "Missing matrix height in Matrix Market file."
              << std::endl;
    return false;
  }

  // Determine the number of columns.
  if (description.object == kMatrixMarketObjectMatrix) {
    if (!(line_stream >> *num_columns)) {
      std::cerr << "Missing matrix width in Matrix Market file."
                << std::endl;
      return false;
    }
  } else {
    *num_columns = 1;
  }

  return true;
}

inline bool ReadMatrixMarketCoordinateMetadata(
    const MatrixMarketDescription& description,
    std::ifstream& file,
    Int* num_rows,
    Int* num_columns,
    Int* num_entries) {
  std::string line;

  // Skip any comment lines.
  while (file.peek() == kMatrixMarketCommentChar) {
    std::getline(file, line);
  }

  // Get a stringstream for the relevant line of the file.
  if (!std::getline(file, line)) {
    std::cerr << "Could not extract the coordinate metadata line." << std::endl;
    return false;
  }
  std::stringstream line_stream(line);

  // Read the number of rows.
  if (!(line_stream >> *num_rows)) {
    std::cerr << "Missing matrix height in Matrix Market file." << std::endl;
    return false;
  }

  // Determine the number of columns.
  if (description.object == kMatrixMarketObjectMatrix) {
    if (!(line_stream >> *num_columns)) {
      std::cerr << "Missing matrix width in Matrix Market file." << std::endl;
      return false;
    }
  } else {
    *num_columns = 1;
  }

  // Read the number of entries.
  if (!(line_stream >> *num_entries)) {
    std::cerr << "Missing num_nonzeros in Matrix Market file." << std::endl;
    return false;
  }

  return true;
}

inline bool ReadMatrixMarketArrayRealValue(
    const MatrixMarketDescription& description,
    std::ifstream& file,
    double* value) {
  std::string line;

  // Skip any comment lines.
  while (file.peek() == kMatrixMarketCommentChar) {
    std::getline(file, line);
  }

  // Get a stringstream for the relevant line of the file.
  if (!std::getline(file, line)) {
    std::cerr << "Could not extract dense entry." << std::endl;
    return false;
  }
  std::stringstream line_stream(line);

  // Read the value.
  if (!(line_stream >> *value)) {
    std::cerr << "Could not extract dense entry." << std::endl;
    return false;
  }

  return true;
}

inline bool ReadMatrixMarketArrayComplexValue(
    const MatrixMarketDescription& description,
    std::ifstream& file,
    double* real_value,
    double* imag_value) {
  std::string line;

  // Skip any comment lines.
  while (file.peek() == kMatrixMarketCommentChar) {
    std::getline(file, line);
  }

  // Get a stringstream for the relevant line of the file.
  if (!std::getline(file, line)) {
    std::cerr << "Could not extract dense entry." << std::endl;
    return false;
  }
  std::stringstream line_stream(line);

  // Read the real value.
  if (!(line_stream >> *real_value)) {
    std::cerr << "Could not extract dense entry." << std::endl;
    return false;
  }

  // Read the imaginary value.
  if (!(line_stream >> *imag_value)) {
    std::cerr << "Could not extract dense entry." << std::endl;
    return false;
  }

  return true;
}

inline bool ReadMatrixMarketCoordinateIndices(
    const MatrixMarketDescription& description,
    std::ifstream& file,
    Int* row,
    Int* column) {
  std::string line;

  // Skip any comment lines.
  while (file.peek() == kMatrixMarketCommentChar) {
    std::getline(file, line);
  }

  // Get a stringstream for the relevant line.
  if (!std::getline(file, line)) {
    std::cerr << "Could not extract entry description from Matrix Market "
                 "file." << std::endl;
    return false;
  }
  std::stringstream line_stream(line);

  // Read the row index.
  if (!(line_stream >> *row)) {
    std::cerr << "Could not extract row index of entry." << std::endl;
    return false;
  }
  --(*row); // Convert from 1-based to 0-based indexing.

  // Determine the column index.
  if (description.object == kMatrixMarketObjectMatrix) {
    if (!(line_stream >> *column)) {
      std::cerr << "Could not extract column index of entry." << std::endl;
      return false;
    }
    --(*column); // Convert from 1-based to 0-based indexing.
  } else {
    *column = 0;
  }

  return true;
}

inline bool ReadMatrixMarketCoordinateRealEntry(
    const MatrixMarketDescription& description,
    std::ifstream& file,
    Int* row,
    Int* column,
    double* value) {
  std::string line;

  // Skip any comment lines.
  while (file.peek() == kMatrixMarketCommentChar) {
    std::getline(file, line);
  }

  // Get a stringstream for the relevant line.
  if (!std::getline(file, line)) {
    std::cerr << "Could not extract entry description from Matrix Market "
                 "file." << std::endl;
    return false;
  }
  std::stringstream line_stream(line);

  // Read the row index.
  if (!(line_stream >> *row)) {
    std::cerr << "Could not extract row index of entry." << std::endl;
    return false;
  }
  --(*row); // Convert from 1-based to 0-based indexing.

  // Determine the column index.
  if (description.object == kMatrixMarketObjectMatrix) {
    if (!(line_stream >> *column)) {
      std::cerr << "Could not extract column index of entry." << std::endl;
      return false;
    }
    --(*column); // Convert from 1-based to 0-based indexing.
  } else {
    *column = 0;
  }

  // Determine the value.
  if (description.field == kMatrixMarketFieldPattern) {
    *value = 1.;
  } else {
    if (!(line_stream >> *value)) {
      std::cerr << "Could not extract value of entry." << std::endl;
      return false;
    }
  }

  return true;
}

inline bool ReadMatrixMarketCoordinateComplexEntry(
    const MatrixMarketDescription& description,
    std::ifstream& file,
    Int* row,
    Int* column,
    double* real_value,
    double* imag_value) {
  std::string line;

  // Skip any comment lines.
  while (file.peek() == kMatrixMarketCommentChar) {
    std::getline(file, line);
  }

  // Get a stringstream for the relevant line.
  if (!std::getline(file, line)) {
    std::cerr << "Could not extract entry description from Matrix Market "
                 "file." << std::endl;
    return false;
  }
  std::stringstream line_stream(line);

  // Read the row index.
  if (!(line_stream >> *row)) {
    std::cerr << "Could not extract row index of entry." << std::endl;
    return false;
  }
  --(*row); // Convert from 1-based to 0-based indexing.

  // Determine the column index.
  if (description.object == kMatrixMarketObjectMatrix) {
    if (!(line_stream >> *column)) {
      std::cerr << "Could not extract column index of entry." << std::endl;
      return false;
    }
    --(*column); // Convert from 1-based to 0-based indexing.
  } else {
    *column = 0;
  }

  // Determine the value.
  if (description.field == kMatrixMarketFieldPattern) {
    *real_value = 1.;
    *imag_value = 0.;
  } else {
    if (!(line_stream >> *real_value)) {
      std::cerr << "Could not extract real value of entry." << std::endl;
      return false;
    }
    if (!(line_stream >> *imag_value)) {
      std::cerr << "Could not extract imag value of entry." << std::endl;
      return false;
    }
  }

  return true;
}

} // namespace quotient

#endif // ifndef QUOTIENT_MATRIX_MARKET_IMPL_H_
