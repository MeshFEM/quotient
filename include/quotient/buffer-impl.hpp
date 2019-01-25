/*
 * Copyright (c) 2019 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_BUFFER_IMPL_H_
#define QUOTIENT_BUFFER_IMPL_H_

#include "quotient/buffer.hpp"

namespace quotient {

template <typename T>
Buffer<T>::Buffer() : size_(0), capacity_(0) {}

template <typename T>
Buffer<T>::Buffer(Int num_elements)
    : size_(num_elements),
      capacity_(num_elements),
      data_(new T[num_elements]) {}

template <typename T>
Buffer<T>::Buffer(Int num_elements, const T& value)
    : size_(num_elements), capacity_(num_elements), data_(new T[num_elements]) {
  // TODO(Jack Poulson): Decide if a variant of 'new' would be faster.
  T* data = Data();
  for (Int index = 0; index < num_elements; ++index) {
    data[index] = value;
  }
}

template <typename T>
Buffer<T>::Buffer(std::initializer_list<T> list)
    : size_(list.size()), capacity_(list.size()), data_(new T[list.size()]) {
  T* data = Data();
  const T* list_data = list.begin();
  for (Int index = 0; index < size_; ++index) {
    data[index] = list_data[index];
  }
}

template <typename T>
Buffer<T>::Buffer(const T* begin, const T* end) {
  size_ = std::distance(begin, end);
  capacity_ = size_;
  data_.reset(new T[size_]);

  T* data = data_.get();
  for (Int index = 0; index < size_; ++index) {
    data[index] = begin[index];
  }
}

template <typename T>
Buffer<T>::Buffer(const Buffer<T>& buffer)
    : size_(buffer.size_), capacity_(buffer.size_), data_(new T[buffer.size_]) {
  T* data = Data();
  const T* buffer_data = buffer.Data();
  for (Int index = 0; index < size_; ++index) {
    data[index] = buffer_data[index];
  }
}

template <typename T>
Buffer<T>::Buffer(const std::vector<T>& vec)
    : size_(vec.size()), capacity_(vec.size()), data_(new T[vec.size()]) {
  T* data = Data();
  for (Int index = 0; index < size_; ++index) {
    data[index] = vec[index];
  }
}

template <typename T>
Buffer<T>& Buffer<T>::operator=(const Buffer<T>& buffer) {
  if (this != &buffer) {
    Resize(buffer.Size());
    T* data = Data();
    const T* buffer_data = buffer.Data();
    for (Int index = 0; index < size_; ++index) {
      data[index] = buffer_data[index];
    }
  }
  return *this;
}

template <typename T>
Buffer<T>& Buffer<T>::operator=(const std::vector<T>& vec) {
  Resize(vec.size());
  T* data = Data();
  for (Int index = 0; index < size_; ++index) {
    data[index] = vec[index];
  }
  return *this;
}

template <typename T>
bool Buffer<T>::operator==(const Buffer<T>& buffer) const {
  if (size_ != buffer.Size()) {
    return false;
  }

  const T* data = data_.get();
  for (Int index = 0; index < size_; ++index) {
    if (data[index] != buffer[index]) {
      return false;
    }
  }

  return true;
}

template <typename T>
bool Buffer<T>::operator!=(const Buffer<T>& buffer) const {
  return !this->operator==(buffer);
}

template <typename T>
void Buffer<T>::Resize(Int num_elements) {
  if (num_elements > capacity_) {
    data_.reset(new T[num_elements]);
    capacity_ = num_elements;
  }
  size_ = num_elements;
}

template <typename T>
void Buffer<T>::Resize(Int num_elements, const T& value) {
  if (num_elements > capacity_) {
    data_.reset(new T[num_elements]);
    capacity_ = num_elements;
  }
  size_ = num_elements;

  T* data = data_.get();
  for (Int index = 0; index < num_elements; ++index) {
    data[index] = value;
  }
}

template <typename T>
Int Buffer<T>::Size() const {
  return size_;
}

template <typename T>
bool Buffer<T>::Empty() const {
  return size_ == 0;
}

template <typename T>
Int Buffer<T>::Capacity() const {
  return capacity_;
}

template <typename T>
T* Buffer<T>::Data() {
  return data_.get();
}

template <typename T>
const T* Buffer<T>::Data() const {
  return data_.get();
}

template <typename T>
T* Buffer<T>::begin() {
  return data_.get();
}

template <typename T>
const T* Buffer<T>::begin() const {
  return data_.get();
}

template <typename T>
T* Buffer<T>::end() {
  return data_.get() + size_;
}

template <typename T>
const T* Buffer<T>::end() const {
  return data_.get() + size_;
}

template <typename T>
T& Buffer<T>::operator[](Int index) {
  T* data = data_.get();
  return data[index];
}

template <typename T>
const T& Buffer<T>::operator[](Int index) const {
  const T* data = data_.get();
  return data[index];
}

template <typename T>
void Buffer<T>::Clear() {
  size_ = 0;
  capacity_ = 0;
  data_.reset();
}

}  // namespace quotient

#endif  // ifndef QUOTIENT_BUFFER_IMPL_H_
