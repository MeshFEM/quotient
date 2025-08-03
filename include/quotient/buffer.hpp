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

template<typename T>
struct UnderlyingDataTraits {
    using type = T;
    using size_type = UInt;
    static constexpr size_type underlying_data_size(size_type element_size) { return element_size; }
    static constexpr T        *underlying_data_ptr(T *base) { return base; }

    template<class allocator>
    static void value_construct(allocator &a, T *ptr, const T &value) {
        std::allocator_traits<allocator>::construct(a, ptr, value);
    }

    template<class allocator>
    static T *allocate(allocator &a, size_type num_elements) {
        return std::allocator_traits<allocator>::allocate(a, num_elements);
    }
};

// This specialization exists because std::complex does *not* satisfy
// std::is_trivially_constructible or std::is_trivially_copy_constructible,
// but 'float' and 'double' do.
// We store the data of `Buffer<Complex<Real>>` as an array of `Real` with
// twice as many elements.
template<typename Real>
struct UnderlyingDataTraits<Complex<Real>> {
    using type = Real;
    using size_type = UInt;
    static_assert(sizeof(Complex<Real>) == 2 * sizeof(Real), "Complex data type should consist of two underlying reals");
    static constexpr size_type underlying_data_size (size_type element_size) { return 2 * element_size; }
    static constexpr type     *underlying_data_ptr  (Complex<Real> *base) { return reinterpret_cast<type *>(base); }

    template<class allocator>
    static void value_construct(allocator &a, Complex<Real> *ptr, const Complex<Real> &value) {
        std::allocator_traits<allocator>::construct(a, underlying_data_ptr(ptr)    , value.real());
        std::allocator_traits<allocator>::construct(a, underlying_data_ptr(ptr) + 1, value.imag());
    }

    template<class allocator>
    static auto allocate(allocator &a, size_type num_elements) {
        return reinterpret_cast<Complex<Real> *>(std::allocator_traits<allocator>::allocate(a, underlying_data_size(num_elements)));
    }
};

// A simple stand-in for some of the functionality of std::vector that avoids
// explicit initialization in the single-argument versions of the constructor
// and in resize. In exchange, it does not have 'push_back' functionality.
template <typename T, class Allocator_ = std::allocator<typename UnderlyingDataTraits<T>::type>>
class Buffer {
 public:
  typedef UnderlyingDataTraits<T> UDT;

  typedef Allocator_ Allocator;
  typedef std::allocator_traits<Allocator> AllocatorTraits;
  typedef typename UDT::size_type SizeType;
  typedef typename UDT::type UnderlyingType;
  typedef T* Pointer;
  typedef const T* ConstPointer;
  typedef Pointer Iterator;
  typedef ConstPointer ConstIterator;
  typedef T& Reference;
  typedef const T& ConstReference;
  static constexpr bool is_trivially_constructible      = std::is_trivially_constructible_v<typename UDT::type>;
  static constexpr bool is_trivially_copy_constructible = std::is_trivially_copy_constructible_v<typename UDT::type>;
  static constexpr bool is_trivially_destructible       = std::is_trivially_destructible_v<typename UDT::type>;

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
  Buffer(const Buffer& buffer);

  // A move constructor.
  Buffer(Buffer&& buffer) QUOTIENT_NOEXCEPT;

  // A copy constructor from std::vector<T>.
  Buffer(const std::vector<T>& vec);

  // An assignment operator.
  Buffer& operator=(const Buffer& buffer);

  // A move assignment operator.
  Buffer& operator=(Buffer&& buffer) QUOTIENT_NOEXCEPT;

  // An assignment operator from std::vector<T>.
  Buffer& operator=(const std::vector<T>& buffer);

  // An equality check against another buffer.
  bool operator==(const Buffer& buffer) const;

  // An inequality check against another buffer.
  bool operator!=(const Buffer& buffer) const;

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

#include "quotient/buffer-impl.hpp"

#endif  // ifndef QUOTIENT_BUFFER_H_
