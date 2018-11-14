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
#include <iostream>
#include <sstream>

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

} // namespace quotient

#endif // ifndef QUOTIENT_MATRIX_MARKET_IMPL_H_
