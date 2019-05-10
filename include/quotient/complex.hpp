/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_COMPLEX_H_
#define QUOTIENT_COMPLEX_H_

#include <complex>
#include <type_traits>

#include "quotient/enable_if.hpp"
#include "quotient/macros.hpp"

namespace quotient {

// An extension of std::complex beyond float and double.
template <class Real>
class Complex {
 public:
  // The real and imaginary components of the complex variable.
  Real real, imag;
};

// A specialization of Complex to an underlying real type of 'float'.
template <>
class Complex<float> : public std::complex<float> {
 public:
  // The underlying real type of the complex class.
  typedef float RealType;

  // Imports of the std::complex operators.
  using std::complex<RealType>::operator=;
  using std::complex<RealType>::operator-=;
  using std::complex<RealType>::operator+=;
  using std::complex<RealType>::operator*=;
  using std::complex<RealType>::operator/=;

  // The default constructor.
  inline Complex() QUOTIENT_NOEXCEPT;

  // A copy constructor from a Complex<Real> variable.
  inline Complex(const Complex<RealType>& input) QUOTIENT_NOEXCEPT;

  // A copy constructor from a std::complex variable.
  inline Complex(const std::complex<RealType>& input) QUOTIENT_NOEXCEPT;

  // A copy constructor from a real variable.
  template <class RealInputType>
  Complex(const RealInputType& input) QUOTIENT_NOEXCEPT;

  // A copy constructor from real and imaginary parts.
  template <class RealInputType, class ImagInputType>
  Complex(const RealInputType& real,
          const ImagInputType& imag) QUOTIENT_NOEXCEPT;

  // A copy constructor from a Complex variable.
  template <class RealInputType>
  Complex(const Complex<RealInputType>& input) QUOTIENT_NOEXCEPT;
};

// A specialization of Complex to an underlying real type of 'double'.
template <>
class Complex<double> : public std::complex<double> {
 public:
  // The underlying real type of the complex class.
  typedef double RealType;

  // Imports of the std::complex operators.
  using std::complex<RealType>::operator=;
  using std::complex<RealType>::operator-=;
  using std::complex<RealType>::operator+=;
  using std::complex<RealType>::operator*=;
  using std::complex<RealType>::operator/=;

  // The default constructor.
  inline Complex() QUOTIENT_NOEXCEPT;

  // A copy constructor from a Complex<Real> variable.
  inline Complex(const Complex<RealType>& input) QUOTIENT_NOEXCEPT;

  // A copy constructor from a std::complex variable.
  inline Complex(const std::complex<RealType>& input) QUOTIENT_NOEXCEPT;

  // A copy constructor from a real variable.
  template <class RealInputType>
  Complex(const RealInputType& input) QUOTIENT_NOEXCEPT;

  // A copy constructor from real and imaginary parts.
  template <class RealInputType, class ImagInputType>
  Complex(const RealInputType& real,
          const ImagInputType& imag) QUOTIENT_NOEXCEPT;

  // A copy constructor from a Complex variable.
  template <class RealInputType>
  Complex(const Complex<RealInputType>& input) QUOTIENT_NOEXCEPT;
};

namespace complex_base {

template <class Real>
struct ComplexBaseHelper {
  typedef Real type;
};

template <class Real>
struct ComplexBaseHelper<Complex<Real>> {
  typedef Real type;
};

}  // namespace complex_base

// Returns the type of the base field of a real or complex scalar. For example:
//   ComplexBase<double> == double
//   ComplexBase<Complex<double>> == double.
template <class Field>
using ComplexBase = typename complex_base::ComplexBaseHelper<Field>::type;

// Encodes whether or not a given type is complex. For example,
//   IsComplex<double>::value == false
//   IsComplex<Complex<double>>::value == true
template <class Real>
struct IsComplex {
  static constexpr bool value = false;
};

template <class Real>
struct IsComplex<Complex<Real>> {
  static constexpr bool value = true;
};

// Encodes whether or not a given type is real. For example,
//   IsComplex<double>::value == true
//   IsComplex<Complex<double>>::value == false
template <class Field>
struct IsReal {
  static constexpr bool value = !IsComplex<Field>::value;
};

// Returns the negation of a complex value.
template <class Real>
Complex<Real> operator-(const Complex<Real>& value) QUOTIENT_NOEXCEPT;

// Returns the sum of two values.
template <class Real>
Complex<Real> operator+(const Complex<Real>& a,
                        const Complex<Real>& b) QUOTIENT_NOEXCEPT;
template <class Real>
Complex<Real> operator+(const Complex<Real>& a,
                        const Real& b) QUOTIENT_NOEXCEPT;
template <class Real>
Complex<Real> operator+(const Real& a,
                        const Complex<Real>& b) QUOTIENT_NOEXCEPT;

// Returns the difference of two values.
template <class Real>
Complex<Real> operator-(const Complex<Real>& a,
                        const Complex<Real>& b) QUOTIENT_NOEXCEPT;
template <class Real>
Complex<Real> operator-(const Complex<Real>& a,
                        const Real& b) QUOTIENT_NOEXCEPT;
template <class Real>
Complex<Real> operator-(const Real& a,
                        const Complex<Real>& b) QUOTIENT_NOEXCEPT;

// Returns the product of two values.
template <class Real>
Complex<Real> operator*(const Complex<Real>& a,
                        const Complex<Real>& b)QUOTIENT_NOEXCEPT;
template <class Real>
Complex<Real> operator*(const Complex<Real>& a, const Real& b)QUOTIENT_NOEXCEPT;
template <class Real>
Complex<Real> operator*(const Real& a, const Complex<Real>& b)QUOTIENT_NOEXCEPT;

// Returns the ratio of two values.
template <class Real>
Complex<Real> operator/(const Complex<Real>& a,
                        const Complex<Real>& b) QUOTIENT_NOEXCEPT;
template <class Real>
Complex<Real> operator/(const Complex<Real>& a,
                        const Real& b) QUOTIENT_NOEXCEPT;
template <class Real>
Complex<Real> operator/(const Real& a,
                        const Complex<Real>& b) QUOTIENT_NOEXCEPT;

// Returns the real part of a real scalar.
template <class Real>
Real RealPart(const Real& value) QUOTIENT_NOEXCEPT;

// Returns the real part of a complex scalar.
template <class Real>
Real RealPart(const Complex<Real>& value) QUOTIENT_NOEXCEPT;

// Returns the imaginary part of a real scalar (zero).
template <class Real>
Real ImagPart(const Real& value) QUOTIENT_NOEXCEPT;

// Returns the imaginary part of a complex scalar.
template <class Real>
Real ImagPart(const Complex<Real>& value) QUOTIENT_NOEXCEPT;

// Returns the complex-conjugate of a real value (the value itself).
template <class Real>
Real Conjugate(const Real& value) QUOTIENT_NOEXCEPT;

// Returns the complex-conjugate of a complex value.
template <class Real>
Complex<Real> Conjugate(const Complex<Real>& value) QUOTIENT_NOEXCEPT;

}  // namespace quotient

#include "quotient/complex-impl.hpp"

#endif  // ifndef QUOTIENT_COMPLEX_H_
