/*
 * Copyright (c) 2019 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_BUFFER_IMPL_H_
#define QUOTIENT_BUFFER_IMPL_H_

#include <algorithm>

#include "quotient/buffer.hpp"

namespace quotient {

template <typename T>
inline Buffer<T>::Buffer() QUOTIENT_NOEXCEPT : size_(0),
                                               capacity_(0),
                                               data_(nullptr) {}

template <typename T>
inline void Buffer<T>::DestructData() {
  DestructRange(0, size_);
  AllocatorTraits::deallocate(allocator_, data_, capacity_);
}

template <typename T>
inline void Buffer<T>::DestructRange(SizeType start, SizeType end) {
  if (!is_trivially_destructible) {
    for (Iterator iter = data_ + start; iter != data_ + end; ++iter) {
      iter->~T();
    }
  }
}

template <typename T>
inline Buffer<T>::~Buffer() {
  DestructData();
}

template <typename T>
inline void Buffer<T>::ConstructRange(SizeType start, SizeType end) {
  if (!is_trivially_constructible) {
    for (Iterator iter = data_ + start; iter != data_ + end; ++iter) {
      AllocatorTraits::construct(allocator_, iter);
    }
  }
}

template <typename T>
inline Buffer<T>::Buffer(SizeType num_elements)
    : size_(num_elements), capacity_(num_elements) {
  data_ = AllocatorTraits::allocate(allocator_, capacity_);
  ConstructRange(0, size_);
}

template <typename T>
inline void Buffer<T>::FillConstructRange(SizeType start, SizeType end,
                                          ConstReference value) {
  if (is_trivially_copy_constructible) {
    std::fill(data_ + start, data_ + end, value);
  } else {
    // Contruct the active elements in-place.
    for (Pointer iter = data_ + start; iter != data_ + end; ++iter) {
      AllocatorTraits::construct(allocator_, iter, value);
    }
  }
}

template <typename T>
inline void Buffer<T>::CopyConstructRange(SizeType start, ConstIterator begin,
                                          ConstIterator end) {
  if (is_trivially_copy_constructible) {
    std::copy(begin, end, data_ + start);
  } else {
    SizeType offset = start;
    for (ConstIterator iter = begin; iter != end; ++iter, ++offset) {
      AllocatorTraits::construct(allocator_, data_ + offset, *iter);
    }
  }
}

template <typename T>
inline Buffer<T>::Buffer(SizeType num_elements, ConstReference value)
    : size_(num_elements), capacity_(num_elements) {
  data_ = AllocatorTraits::allocate(allocator_, capacity_);
  FillConstructRange(0, size_, value);
}

template <typename T>
Buffer<T>::Buffer(ConstIterator begin, ConstIterator end) {
  size_ = std::distance(begin, end);
  capacity_ = size_;
  data_ = AllocatorTraits::allocate(allocator_, capacity_);
  CopyConstructRange(0, begin, end);
}

template <typename T>
Buffer<T>::Buffer(std::initializer_list<T> list)
    : Buffer(list.begin(), list.end()) {}

template <typename T>
inline Buffer<T>::Buffer(const Buffer<T>& buffer)
    : Buffer(buffer.begin(), buffer.end()) {}

template <typename T>
inline Buffer<T>::Buffer(const std::vector<T>& vec)
    : Buffer(vec.begin(), vec.end()) {}

template <typename T>
inline Buffer<T>::Buffer(Buffer<T>&& buffer) QUOTIENT_NOEXCEPT
    : size_(buffer.size_),
      capacity_(buffer.capacity_),
      data_(buffer.data_) {
  buffer.size_ = 0;
  buffer.capacity_ = 0;
  buffer.data_ = nullptr;
}

template <typename T>
Buffer<T>& Buffer<T>::operator=(const Buffer<T>& buffer) {
  if (this != &buffer) {
    const SizeType num_elements = buffer.Size();
    if (num_elements > capacity_) {
      DestructData();
      data_ = AllocatorTraits::allocate(allocator_, num_elements);
      size_ = num_elements;
      capacity_ = num_elements;
      CopyConstructRange(0, buffer.begin(), buffer.end());
    } else if (num_elements >= size_) {
      std::copy(buffer.begin(), buffer.begin() + size_, data_);
      CopyConstructRange(size_, buffer.begin() + size_, buffer.end());
      size_ = num_elements;
    } else {
      DestructRange(num_elements, size_);
      size_ = num_elements;
      std::copy(buffer.begin(), buffer.end(), data_);
    }
  }
  return *this;
}

template <typename T>
Buffer<T>& Buffer<T>::operator=(Buffer<T>&& buffer) QUOTIENT_NOEXCEPT {
  DestructData();

  size_ = buffer.size_;
  capacity_ = buffer.capacity_;
  data_ = buffer.data_;

  buffer.size_ = 0;
  buffer.capacity_ = 0;
  buffer.data_ = nullptr;

  return *this;
}

template <typename T>
Buffer<T>& Buffer<T>::operator=(const std::vector<T>& vec) {
  const SizeType num_elements = vec.size();
  if (num_elements > capacity_) {
    DestructData();
    data_ = AllocatorTraits::allocate(allocator_, num_elements);
    size_ = num_elements;
    capacity_ = num_elements;
    CopyConstructRange(0, vec.data(), vec.data() + num_elements);
  } else if (num_elements >= size_) {
    std::copy(vec.begin(), vec.begin() + size_, data_);
    CopyConstructRange(size_, vec.data() + size_, vec.data() + num_elements);
    size_ = num_elements;
  } else {
    DestructRange(num_elements, size_);
    size_ = num_elements;
    std::copy(vec.begin(), vec.end(), data_);
  }

  return *this;
}

template <typename T>
bool Buffer<T>::operator==(const Buffer<T>& buffer) const {
  if (size_ != buffer.Size()) {
    return false;
  }

  for (SizeType index = 0; index < size_; ++index) {
    if (data_[index] != buffer[index]) {
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
void Buffer<T>::Resize(SizeType num_elements) {
  if (num_elements > capacity_) {
    DestructData();
    data_ = AllocatorTraits::allocate(allocator_, num_elements);
    size_ = num_elements;
    capacity_ = num_elements;
    ConstructRange(0, size_);
  } else if (num_elements >= size_) {
    ConstructRange(size_, num_elements);
    size_ = num_elements;
  } else {
    DestructRange(num_elements, size_);
    size_ = num_elements;
  }
}

template <typename T>
void Buffer<T>::Resize(SizeType num_elements, ConstReference value) {
  if (num_elements > capacity_) {
    DestructData();
    data_ = AllocatorTraits::allocate(allocator_, num_elements);
    size_ = num_elements;
    capacity_ = num_elements;
    FillConstructRange(0, size_, value);
  } else if (num_elements >= size_) {
    std::fill(data_, data_ + size_, value);
    FillConstructRange(size_, num_elements, value);
    size_ = num_elements;
  } else {
    DestructRange(num_elements, size_);
    size_ = num_elements;
    std::fill(data_, data_ + size_, value);
  }
}

template <typename T>
inline typename Buffer<T>::SizeType Buffer<T>::Size() const QUOTIENT_NOEXCEPT {
  return size_;
}

template <typename T>
inline bool Buffer<T>::Empty() const QUOTIENT_NOEXCEPT {
  return size_ == 0;
}

template <typename T>
inline typename Buffer<T>::SizeType Buffer<T>::Capacity() const
    QUOTIENT_NOEXCEPT {
  return capacity_;
}

template <typename T>
inline typename Buffer<T>::Pointer Buffer<T>::Data() QUOTIENT_NOEXCEPT {
  return data_;
}

template <typename T>
inline typename Buffer<T>::ConstPointer Buffer<T>::Data() const
    QUOTIENT_NOEXCEPT {
  return data_;
}

template <typename T>
inline typename Buffer<T>::Pointer Buffer<T>::begin() QUOTIENT_NOEXCEPT {
  return data_;
}

template <typename T>
inline typename Buffer<T>::ConstPointer Buffer<T>::begin() const
    QUOTIENT_NOEXCEPT {
  return data_;
}

template <typename T>
inline typename Buffer<T>::ConstPointer Buffer<T>::cbegin() const
    QUOTIENT_NOEXCEPT {
  return data_;
}

template <typename T>
inline typename Buffer<T>::Pointer Buffer<T>::end() QUOTIENT_NOEXCEPT {
  return data_ + size_;
}

template <typename T>
inline typename Buffer<T>::ConstPointer Buffer<T>::end() const
    QUOTIENT_NOEXCEPT {
  return data_ + size_;
}

template <typename T>
inline typename Buffer<T>::ConstPointer Buffer<T>::cend() const
    QUOTIENT_NOEXCEPT {
  return data_ + size_;
}

template <typename T>
inline typename Buffer<T>::Reference Buffer<T>::operator[](SizeType index)
    QUOTIENT_NOEXCEPT {
  return data_[index];
}

template <typename T>
inline typename Buffer<T>::ConstReference Buffer<T>::operator[](
    SizeType index) const QUOTIENT_NOEXCEPT {
  return data_[index];
}

template <typename T>
inline typename Buffer<T>::Reference Buffer<T>::Back() QUOTIENT_NOEXCEPT {
  return data_[size_ - 1];
}

template <typename T>
inline typename Buffer<T>::ConstReference Buffer<T>::Back() const
    QUOTIENT_NOEXCEPT {
  return data_[size_ - 1];
}

template <typename T>
inline void Buffer<T>::Clear() QUOTIENT_NOEXCEPT {
  DestructData();
  size_ = 0;
  capacity_ = 0;
  data_ = nullptr;
}

}  // namespace quotient

#endif  // ifndef QUOTIENT_BUFFER_IMPL_H_
