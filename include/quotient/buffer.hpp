/*
 * Copyright (c) 2019 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_BUFFER_H_
#define QUOTIENT_BUFFER_H_

#include <initializer_list>
#include <memory>
#include <vector>

#include "quotient/complex.hpp"
#include "quotient/integers.hpp"
#include "quotient/macros.hpp"

namespace quotient {

// A simple stand-in for some of the functionality of std::vector that avoids
// explicit initialization in the single-argument versions of the constructor
// and in resize. In exchange, it does not have 'push_back' functionality.
template <typename T>
class Buffer {
 public:
  typedef std::allocator<T> Allocator;
  typedef std::allocator_traits<Allocator> AllocatorTraits;
  typedef UInt SizeType;
  typedef T* Pointer;
  typedef const T* ConstPointer;
  typedef Pointer Iterator;
  typedef ConstPointer ConstIterator;
  typedef T& Reference;
  typedef const T& ConstReference;
  static constexpr bool is_trivially_constructible =
      std::is_trivially_constructible<T>::value;
  static constexpr bool is_trivially_copy_constructible =
      std::is_trivially_copy_constructible<T>::value;
  static constexpr bool is_trivially_destructible =
      std::is_trivially_destructible<T>::value;

  // Constructs a zero-length buffer.
  Buffer() QUOTIENT_NOEXCEPT;

  // Destructor.
  ~Buffer();

  // Constructs a buffer of the given length without initializing the data.
  Buffer(SizeType num_elements);

  // Constructs a buffer of the given length where each entry is initialized
  // to the specified value.
  Buffer(SizeType num_elements, ConstReference value);

  // Constructs the buffer from an initializer list.
  Buffer(std::initializer_list<T> list);

  // Constructs a buffer by copying the given iterator range.
  Buffer(ConstPointer begin, ConstPointer end);

  // A copy constructor.
  Buffer(const Buffer<T>& buffer);

  // A move constructor.
  Buffer(Buffer<T>&& buffer) QUOTIENT_NOEXCEPT;

  // A copy constructor from std::vector<T>.
  Buffer(const std::vector<T>& vec);

  // An assignment operator.
  Buffer<T>& operator=(const Buffer<T>& buffer);

  // A move assignment operator.
  Buffer<T>& operator=(Buffer<T>&& buffer) QUOTIENT_NOEXCEPT;

  // An assignment operator from std::vector<T>.
  Buffer<T>& operator=(const std::vector<T>& buffer);

  // An equality check against another buffer.
  bool operator==(const Buffer<T>& buffer) const;

  // An inequality check against another buffer.
  bool operator!=(const Buffer<T>& buffer) const;

  // Resizes the buffer to store the given number of elements and avoids
  // initializing values where possible.
  void Resize(SizeType num_elements);

  // Resizes the buffer to store the given number of elements and initializes
  // *all* entries in the new buffer to the specified value.
  void Resize(SizeType num_elements, ConstReference value);

  // Returns the current length of the buffer.
  SizeType Size() const QUOTIENT_NOEXCEPT;

  // Returns true if the buffer does not contain any items.
  bool Empty() const QUOTIENT_NOEXCEPT;

  // Returns the current capacity of the buffer.
  SizeType Capacity() const QUOTIENT_NOEXCEPT;

  // Returns a mutable pointer to the underlying buffer of entries.
  Pointer Data() QUOTIENT_NOEXCEPT;

  // Returns an immutable pointer to the underlying buffer of entries.
  ConstPointer Data() const QUOTIENT_NOEXCEPT;

  // A mutable iterator to the beginning of the buffer.
  // NOLINTNEXTLINE(readability-identifier-naming)
  Iterator begin() QUOTIENT_NOEXCEPT;

  // An immutable iterator to the beginning of the buffer.
  // NOLINTNEXTLINE(readability-identifier-naming)
  ConstIterator begin() const QUOTIENT_NOEXCEPT;

  // An immutable iterator to the beginning of the buffer.
  // NOLINTNEXTLINE(readability-identifier-naming)
  ConstIterator cbegin() const QUOTIENT_NOEXCEPT;

  // A mutable iterator for the end of the buffer.
  // NOLINTNEXTLINE(readability-identifier-naming)
  Iterator end() QUOTIENT_NOEXCEPT;

  // An immutable iterator for the end of the buffer.
  // NOLINTNEXTLINE(readability-identifier-naming)
  ConstIterator end() const QUOTIENT_NOEXCEPT;

  // An immutable iterator for the end of the buffer.
  // NOLINTNEXTLINE(readability-identifier-naming)
  ConstIterator cend() const QUOTIENT_NOEXCEPT;

  // Returns a mutable reference to the given entry of the buffer.
  Reference operator[](SizeType index) QUOTIENT_NOEXCEPT;

  // Returns an immutable reference to the given entry of the buffer.
  ConstReference operator[](SizeType index) const QUOTIENT_NOEXCEPT;

  // Returns a reference to the last element of the list.
  // Its behavior is undefined if the list is empty.
  Reference Back() QUOTIENT_NOEXCEPT;

  // Returns an immmutable reference to the last element of the list.
  // Its behavior is undefined if the list is empty.
  ConstReference Back() const QUOTIENT_NOEXCEPT;

  // Frees all memory and sets the buffer size to zero.
  void Clear() QUOTIENT_NOEXCEPT;

 private:
  // The number of entries stored by the buffer.
  SizeType size_;

  // The number of entries that could have been stored by 'data_'.
  SizeType capacity_;

  // The underlying pointer for storing entries in the buffer.
  Pointer data_;

  // The allocator to be used for managing memory.
  Allocator allocator_;

  // Destruct the constructed members.
  void DestructData();

  // Destructs the entries in the specified range.
  void DestructRange(SizeType start, SizeType end);

  // Default-constructs the elements in the specified range.
  void ConstructRange(SizeType start, SizeType end);

  // Constructs the elements in the specified positional range by
  // copy-constructing with the specified value.
  void FillConstructRange(SizeType start, SizeType end, ConstReference value);

  // Construct the elements starting at the given index by traversing the
  // specified range of inputs.
  void CopyConstructRange(SizeType start, ConstIterator begin,
                          ConstIterator end);
};

// This specialization exists because std::complex does *not* satisfy
// std::is_trivially_constructible or std::is_trivially_copy_constructible,
// but 'float' and 'double' do.
template <typename Real>
class Buffer<Complex<Real>> {
 public:
  typedef std::allocator<Real> Allocator;
  typedef std::allocator_traits<Allocator> AllocatorTraits;
  typedef UInt SizeType;
  typedef Complex<Real>* Pointer;
  typedef const Complex<Real>* ConstPointer;
  typedef Pointer Iterator;
  typedef ConstPointer ConstIterator;
  typedef Complex<Real>& Reference;
  typedef const Complex<Real>& ConstReference;
  static constexpr bool is_trivially_constructible_real =
      std::is_trivially_constructible<Real>::value;
  static constexpr bool is_trivially_copy_constructible_real =
      std::is_trivially_copy_constructible<Real>::value;
  static constexpr bool is_trivially_destructible_real =
      std::is_trivially_destructible<Real>::value;

  // Constructs a zero-length buffer.
  Buffer() QUOTIENT_NOEXCEPT;

  // Destructor.
  ~Buffer();

  // Constructs a buffer of the given length without initializing the data.
  Buffer(SizeType num_elements);

  // Constructs a buffer of the given length where each entry is initialized
  // to the specified value.
  Buffer(SizeType num_elements, ConstReference value);

  // Constructs the buffer from an initializer list.
  Buffer(std::initializer_list<Complex<Real>> list);

  // Constructs a buffer by copying the given iterator range.
  Buffer(ConstPointer begin, ConstPointer end);

  // A copy constructor.
  Buffer(const Buffer<Complex<Real>>& buffer);

  // A move constructor.
  Buffer(Buffer<Complex<Real>>&& buffer) QUOTIENT_NOEXCEPT;

  // A copy constructor from std::vector<Complex<Real>>.
  Buffer(const std::vector<Complex<Real>>& vec);

  // An assignment operator.
  Buffer<Complex<Real>>& operator=(const Buffer<Complex<Real>>& buffer);

  // A move assignment operator.
  Buffer<Complex<Real>>& operator=(Buffer<Complex<Real>>&& buffer)
      QUOTIENT_NOEXCEPT;

  // An assignment operator from std::vector<Complex<Real>>.
  Buffer<Complex<Real>>& operator=(const std::vector<Complex<Real>>& buffer);

  // An equality check against another buffer.
  bool operator==(const Buffer<Complex<Real>>& buffer) const;

  // An inequality check against another buffer.
  bool operator!=(const Buffer<Complex<Real>>& buffer) const;

  // Resizes the buffer to store the given number of elements and avoids
  // initializing values where possible.
  void Resize(SizeType num_elements);

  // Resizes the buffer to store the given number of elements and initializes
  // *all* entries in the new buffer to the specified value.
  void Resize(SizeType num_elements, ConstReference value);

  // Returns the current length of the buffer.
  SizeType Size() const QUOTIENT_NOEXCEPT;

  // Returns true if the buffer does not contain any items.
  bool Empty() const QUOTIENT_NOEXCEPT;

  // Returns the current capacity of the buffer.
  SizeType Capacity() const QUOTIENT_NOEXCEPT;

  // Returns a mutable pointer to the underlying buffer of entries.
  Pointer Data() QUOTIENT_NOEXCEPT;

  // Returns an immutable pointer to the underlying buffer of entries.
  ConstPointer Data() const QUOTIENT_NOEXCEPT;

  // A mutable iterator to the beginning of the buffer.
  // NOLINTNEXTLINE(readability-identifier-naming)
  Iterator begin() QUOTIENT_NOEXCEPT;

  // An immutable iterator to the beginning of the buffer.
  // NOLINTNEXTLINE(readability-identifier-naming)
  ConstIterator begin() const QUOTIENT_NOEXCEPT;

  // An immutable iterator to the beginning of the buffer.
  // NOLINTNEXTLINE(readability-identifier-naming)
  ConstIterator cbegin() const QUOTIENT_NOEXCEPT;

  // A mutable iterator for the end of the buffer.
  // NOLINTNEXTLINE(readability-identifier-naming)
  Iterator end() QUOTIENT_NOEXCEPT;

  // An immutable iterator for the end of the buffer.
  // NOLINTNEXTLINE(readability-identifier-naming)
  ConstIterator end() const QUOTIENT_NOEXCEPT;

  // An immutable iterator for the end of the buffer.
  // NOLINTNEXTLINE(readability-identifier-naming)
  ConstIterator cend() const QUOTIENT_NOEXCEPT;

  // Returns a mutable reference to the given entry of the buffer.
  Reference operator[](SizeType index) QUOTIENT_NOEXCEPT;

  // Returns an immutable reference to the given entry of the buffer.
  ConstReference operator[](SizeType index) const QUOTIENT_NOEXCEPT;

  // Returns a reference to the last element of the list.
  // Its behavior is undefined if the list is empty.
  Reference Back() QUOTIENT_NOEXCEPT;

  // Returns an immmutable reference to the last element of the list.
  // Its behavior is undefined if the list is empty.
  ConstReference Back() const QUOTIENT_NOEXCEPT;

  // Frees all memory and sets the buffer size to zero.
  void Clear() QUOTIENT_NOEXCEPT;

 private:
  // The number of entries stored by the buffer.
  SizeType size_;

  // The number of entries that could have been stored by 'data_'.
  SizeType capacity_;

  // The underlying pointer for storing entries in the buffer.
  Pointer data_;

  // The allocator to be used for managing memory.
  Allocator allocator_;

  // Destruct the constructed members.
  void DestructData();

  // Destructs the entries in the specified range.
  void DestructRange(SizeType start, SizeType end);

  // Default-constructs the elements in the specified range.
  void ConstructRange(SizeType start, SizeType end);

  // Constructs the elements in the specified positional range by
  // copy-constructing with the specified value.
  void FillConstructRange(SizeType start, SizeType end, ConstReference value);

  // Construct the elements starting at the given index by traversing the
  // specified range of inputs.
  void CopyConstructRange(SizeType start, ConstIterator begin,
                          ConstIterator end);
};

}  // namespace quotient

#include "quotient/buffer-complex-impl.hpp"
#include "quotient/buffer-impl.hpp"

#endif  // ifndef QUOTIENT_BUFFER_H_
