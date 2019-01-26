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

#include "quotient/integers.hpp"

namespace quotient {

// A simple stand-in for some of the functionality of std::vector that avoids
// explicit initialization in the single-argument versions of the constructor
// and in resize. In exchange, it does not have 'push_back' functionality.
template <typename T>
class Buffer {
 public:
  // Constructs a zero-length buffer.
  Buffer();

  // Constructs a buffer of the given length without initializing the data.
  Buffer(Int num_elements);

  // Constructs a buffer of the given length where each entry is initialized
  // to the specified value.
  Buffer(Int num_elements, const T& value);

  // Constructs the buffer from an initializer list.
  Buffer(std::initializer_list<T> list);

  // Constructs a buffer by copying the given iterator range.
  Buffer(const T* begin, const T* end);

  // A copy constructor.
  Buffer(const Buffer<T>& buffer);

  // A copy constructor from std::vector<T>.
  Buffer(const std::vector<T>& vec);

  // An assignment operator.
  Buffer<T>& operator=(const Buffer<T>& buffer);

  // An assignment operator from std::vector<T>.
  Buffer<T>& operator=(const std::vector<T>& buffer);

  // An equality check against another buffer.
  bool operator==(const Buffer<T>& buffer) const;

  // An inequality check against another buffer.
  bool operator!=(const Buffer<T>& buffer) const;

  // Resizes the buffer to store the given number of elements and avoids
  // initializing values where possible.
  void Resize(Int num_elements);

  // Resizes the buffer to store the given number of elements and initializes
  // *all* entries in the new buffer to the specified value.
  void Resize(Int num_elements, const T& value);

  // Returns the current length of the buffer.
  Int Size() const;

  // Returns true if the buffer does not contain any items.
  bool Empty() const;

  // Returns the current capacity of the buffer.
  Int Capacity() const;

  // Returns a mutable pointer to the underlying buffer of entries.
  T* Data();

  // Returns an immutable pointer to the underlying buffer of entries.
  const T* Data() const;

  // A mutable iterator to the beginning of the buffer.
  T* begin();  // NOLINT(readability-identifier-naming)

  // An immutable iterator to the beginning of the buffer.
  const T* begin() const;  // NOLINT(readability-identifier-naming)

  // An immutable iterator to the beginning of the buffer.
  const T* cbegin() const;  // NOLINT(readability-identifier-naming)

  // A mutable iterator for the end of the buffer.
  T* end();  // NOLINT(readability-identifier-naming)

  // An immutable iterator for the end of the buffer.
  const T* end() const;  // NOLINT(readability-identifier-naming)

  // An immutable iterator for the end of the buffer.
  const T* cend() const;  // NOLINT(readability-identifier-naming)

  // Returns a mutable reference to the given entry of the buffer.
  T& operator[](Int index);

  // Returns an immutable reference to the given entry of the buffer.
  const T& operator[](Int index) const;

  // Returns a reference to the last element of the list.
  // Its behavior is undefined if the list is empty.
  T& Back();

  // Returns an immmutable reference to the last element of the list.
  // Its behavior is undefined if the list is empty.
  const T& Back() const;

  // Frees all memory and sets the buffer size to zero.
  void Clear();

 private:
  // The number of entries stored by the buffer.
  Int size_;

  // The number of entries that could have been stored by 'data_'.
  Int capacity_;

  // The underlying pointer for storing entries in the buffer.
  std::unique_ptr<T[]> data_;
};

}  // namespace quotient

#include "quotient/buffer-impl.hpp"

#endif  // ifndef QUOTIENT_BUFFER_H_
