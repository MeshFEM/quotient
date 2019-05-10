/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_COMPLEX_IMPL_H_
#define QUOTIENT_COMPLEX_IMPL_H_

#include "quotient/complex.hpp"

namespace quotient {

inline Complex<float>::Complex() QUOTIENT_NOEXCEPT : std::complex<float>() {}

inline Complex<double>::Complex() QUOTIENT_NOEXCEPT : std::complex<double>() {}

inline Complex<float>::Complex(const Complex<float>& input) QUOTIENT_NOEXCEPT
    : std::complex<float>(input.real(), input.imag()) {}

inline Complex<double>::Complex(const Complex<double>& input) QUOTIENT_NOEXCEPT
    : std::complex<double>(input.real(), input.imag()) {}

inline Complex<float>::Complex(const std::complex<float>& input)
    QUOTIENT_NOEXCEPT : std::complex<float>(input) {}

inline Complex<double>::Complex(const std::complex<double>& input)
    QUOTIENT_NOEXCEPT : std::complex<double>(input) {}

template <class RealInputType>
Complex<float>::Complex(const RealInputType& input) QUOTIENT_NOEXCEPT
    : std::complex<float>(static_cast<float>(input)) {}

template <class RealInputType>
Complex<double>::Complex(const RealInputType& input) QUOTIENT_NOEXCEPT
    : std::complex<double>(static_cast<double>(input)) {}

template <class RealInputType>
Complex<float>::Complex(const Complex<RealInputType>& input) QUOTIENT_NOEXCEPT
    : std::complex<float>(static_cast<float>(input.real()),
                          static_cast<float>(input.imag())) {}

template <class RealInputType>
Complex<double>::Complex(const Complex<RealInputType>& input) QUOTIENT_NOEXCEPT
    : std::complex<double>(static_cast<double>(input.real()),
                           static_cast<double>(input.imag())) {}

template <class RealInputType, class ImagInputType>
Complex<float>::Complex(const RealInputType& real,
                        const ImagInputType& imag) QUOTIENT_NOEXCEPT
    : std::complex<float>(static_cast<float>(real), static_cast<float>(imag)) {}

template <class RealInputType, class ImagInputType>
Complex<double>::Complex(const RealInputType& real,
                         const ImagInputType& imag) QUOTIENT_NOEXCEPT
    : std::complex<double>(static_cast<double>(real),
                           static_cast<double>(imag)) {}

template <class Real>
Complex<Real> operator-(const Complex<Real>& value) QUOTIENT_NOEXCEPT {
  const std::complex<Real>& value_std =
      static_cast<const std::complex<Real>&>(value);
  return -value_std;
}

template <class Real>
Complex<Real> operator+(const Complex<Real>& a,
                        const Complex<Real>& b) QUOTIENT_NOEXCEPT {
  const std::complex<Real>& a_std = static_cast<const std::complex<Real>&>(a);
  const std::complex<Real>& b_std = static_cast<const std::complex<Real>&>(b);
  return a_std + b_std;
}

template <class Real>
Complex<Real> operator+(const Complex<Real>& a,
                        const Real& b) QUOTIENT_NOEXCEPT {
  const std::complex<Real>& a_std = static_cast<const std::complex<Real>&>(a);
  return a_std + b;
}

template <class Real>
Complex<Real> operator+(const Real& a,
                        const Complex<Real>& b) QUOTIENT_NOEXCEPT {
  const std::complex<Real>& b_std = static_cast<const std::complex<Real>&>(b);
  return a + b_std;
}

template <class Real>
Complex<Real> operator-(const Complex<Real>& a,
                        const Complex<Real>& b) QUOTIENT_NOEXCEPT {
  const std::complex<Real>& a_std = static_cast<const std::complex<Real>&>(a);
  const std::complex<Real>& b_std = static_cast<const std::complex<Real>&>(b);
  return a_std - b_std;
}

template <class Real>
Complex<Real> operator-(const Complex<Real>& a,
                        const Real& b) QUOTIENT_NOEXCEPT {
  const std::complex<Real>& a_std = static_cast<const std::complex<Real>&>(a);
  return a_std - b;
}

template <class Real>
Complex<Real> operator-(const Real& a,
                        const Complex<Real>& b) QUOTIENT_NOEXCEPT {
  const std::complex<Real>& b_std = static_cast<const std::complex<Real>&>(b);
  return a - b_std;
}

template <class Real>
Complex<Real> operator*(const Complex<Real>& a,
                        const Complex<Real>& b)QUOTIENT_NOEXCEPT {
  const std::complex<Real>& a_std = static_cast<const std::complex<Real>&>(a);
  const std::complex<Real>& b_std = static_cast<const std::complex<Real>&>(b);
  return a_std * b_std;
}

template <class Real>
Complex<Real> operator*(const Complex<Real>& a,
                        const Real& b)QUOTIENT_NOEXCEPT {
  const std::complex<Real>& a_std = static_cast<const std::complex<Real>&>(a);
  return a_std * b;
}

template <class Real>
Complex<Real> operator*(const Real& a,
                        const Complex<Real>& b)QUOTIENT_NOEXCEPT {
  const std::complex<Real>& b_std = static_cast<const std::complex<Real>&>(b);
  return a * b_std;
}

template <class Real>
Complex<Real> operator/(const Complex<Real>& a,
                        const Complex<Real>& b) QUOTIENT_NOEXCEPT {
  const std::complex<Real>& a_std = static_cast<const std::complex<Real>&>(a);
  const std::complex<Real>& b_std = static_cast<const std::complex<Real>&>(b);
  return a_std / b_std;
}

template <class Real>
Complex<Real> operator/(const Complex<Real>& a,
                        const Real& b) QUOTIENT_NOEXCEPT {
  const std::complex<Real>& a_std = static_cast<const std::complex<Real>&>(a);
  return a_std / b;
}

template <class Real>
Complex<Real> operator/(const Real& a,
                        const Complex<Real>& b) QUOTIENT_NOEXCEPT {
  const std::complex<Real>& b_std = static_cast<const std::complex<Real>&>(b);
  return a / b_std;
}

template <class Real>
Real RealPart(const Real& value) QUOTIENT_NOEXCEPT {
  return value;
}

template <class Real>
Real RealPart(const Complex<Real>& value) QUOTIENT_NOEXCEPT {
  return value.real();
}

template <class Real>
Real ImagPart(const Real& value) QUOTIENT_NOEXCEPT {
  return 0;
}

template <class Real>
Real ImagPart(const Complex<Real>& value) QUOTIENT_NOEXCEPT {
  return value.imag();
}

template <class Real>
Real Conjugate(const Real& value) QUOTIENT_NOEXCEPT {
  return value;
}

template <class Real>
Complex<Real> Conjugate(const Complex<Real>& value) QUOTIENT_NOEXCEPT {
  return Complex<Real>{value.real(), -value.imag()};
}

}  // namespace quotient

#endif  // ifndef QUOTIENT_COMPLEX_IMPL_H_
