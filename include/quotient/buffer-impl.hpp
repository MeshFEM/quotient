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

template <typename T, class A>
inline Buffer<T, A>::Buffer() QUOTIENT_NOEXCEPT : size_(0),
                                               capacity_(0),
                                               data_(nullptr) {}

template <typename T, class A>
inline void Buffer<T, A>::DestructData() {
  DestructRange(0, size_);
  AllocatorTraits::deallocate(allocator_, reinterpret_cast<UnderlyingType *>(data_), UDT::underlying_data_size(capacity_));
}

template <typename T, class A>
inline void Buffer<T, A>::DestructRange(SizeType start, SizeType end) {
  if (!is_trivially_destructible) {
    UnderlyingType *data = reinterpret_cast<UnderlyingType *>(data_);
    for (UnderlyingType *iter = UDT::underlying_data_ptr(data_ + start); iter != UDT::underlying_data_ptr(data_ + end); ++iter) {
      iter->~UnderlyingType();
    }
  }
}

template <typename T, class A>
inline Buffer<T, A>::~Buffer() {
  DestructData();
}

template <typename T, class A>
inline void Buffer<T, A>::ConstructRange(SizeType start, SizeType end) {
  if (!is_trivially_constructible) {
    for (UnderlyingType *iter = UDT::underlying_data_ptr(data_ + start); iter != UDT::underlying_data_ptr(data_ + end); ++iter) {
      AllocatorTraits::construct(allocator_, iter);
    }
  }
}

template <typename T, class A>
inline Buffer<T, A>::Buffer(SizeType num_elements)
    : size_(num_elements), capacity_(num_elements) {
  data_ = UDT::allocate(allocator_, num_elements);
  ConstructRange(0, size_);
}

template <typename T, class A>
inline void Buffer<T, A>::FillConstructRange(SizeType start, SizeType end,
                                          ConstReference value) {
  if (is_trivially_copy_constructible) {
    std::fill(data_ + start, data_ + end, value);
  } else {
    // Contruct the active elements in-place.
    for (Pointer iter = data_ + start; iter != data_ + end; ++iter) {
      UDT::value_construct(allocator_, iter, value);
    }
  }
}

template <typename T, class A>
inline void Buffer<T, A>::CopyConstructRange(SizeType start, ConstIterator begin,
                                          ConstIterator end) {
  if (is_trivially_copy_constructible) {
    std::copy(begin, end, data_ + start);
  } else {
    SizeType offset = start;
    for (ConstIterator iter = begin; iter != end; ++iter, ++offset) {
        UDT::value_construct(allocator_, data_ + offset, *iter);
    }
  }
}

template <typename T, class A>
inline Buffer<T, A>::Buffer(SizeType num_elements, ConstReference value)
    : size_(num_elements), capacity_(num_elements) {
  data_ = UDT::allocate(allocator_, num_elements);
  FillConstructRange(0, size_, value);
}

template <typename T, class A>
Buffer<T, A>::Buffer(ConstIterator begin, ConstIterator end) {
  size_ = std::distance(begin, end);
  capacity_ = size_;
  data_ = UDT::allocate(allocator_, capacity_);
  CopyConstructRange(0, begin, end);
}

template <typename T, class A>
Buffer<T, A>::Buffer(std::initializer_list<T> list)
    : Buffer(list.begin(), list.end()) {}

template <typename T, class A>
inline Buffer<T, A>::Buffer(const Buffer& buffer)
    : Buffer(buffer.begin(), buffer.end()) {}

template <typename T, class A>
inline Buffer<T, A>::Buffer(const std::vector<T>& vec)
    : Buffer(vec.begin(), vec.end()) {}

template <typename T, class A>
inline Buffer<T, A>::Buffer(Buffer&& buffer) QUOTIENT_NOEXCEPT
    : size_(buffer.size_),
      capacity_(buffer.capacity_),
      data_(buffer.data_) {
  buffer.size_ = 0;
  buffer.capacity_ = 0;
  buffer.data_ = nullptr;
}

template <typename T, class A>
Buffer<T, A>& Buffer<T, A>::operator=(const Buffer& buffer) {
  if (this != &buffer) {
    const SizeType num_elements = buffer.Size();
    if (num_elements > capacity_) {
      DestructData();
      data_ = UDT::allocate(allocator_, num_elements);
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

template <typename T, class A>
Buffer<T, A>& Buffer<T, A>::operator=(Buffer&& buffer) QUOTIENT_NOEXCEPT {
  DestructData();

  size_ = buffer.size_;
  capacity_ = buffer.capacity_;
  data_ = buffer.data_;

  buffer.size_ = 0;
  buffer.capacity_ = 0;
  buffer.data_ = nullptr;

  return *this;
}

template <typename T, class A>
Buffer<T, A>& Buffer<T, A>::operator=(const std::vector<T>& vec) {
  const SizeType num_elements = vec.size();
  if (num_elements > capacity_) {
    DestructData();
    data_ = UDT::allocate(allocator_, num_elements);
    size_ = num_elements;
    capacity_ = num_elements;
    CopyConstructRange(0, vec.begin(), vec.end());
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

template <typename T, class A>
bool Buffer<T, A>::operator==(const Buffer& buffer) const {
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

template <typename T, class A>
bool Buffer<T, A>::operator!=(const Buffer& buffer) const {
  return !this->operator==(buffer);
}

template <typename T, class A>
void Buffer<T, A>::Resize(SizeType num_elements) {
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

template <typename T, class A>
void Buffer<T, A>::Resize(SizeType num_elements, ConstReference value) {
  if (num_elements > capacity_) {
    DestructData();
    data_ = UDT::allocate(allocator_, num_elements);
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

template <typename T, class A>
inline typename Buffer<T, A>::SizeType Buffer<T, A>::Size() const QUOTIENT_NOEXCEPT {
  return size_;
}

template <typename T, class A>
inline bool Buffer<T, A>::Empty() const QUOTIENT_NOEXCEPT {
  return size_ == 0;
}

template <typename T, class A>
inline typename Buffer<T, A>::SizeType Buffer<T, A>::Capacity() const
    QUOTIENT_NOEXCEPT {
  return capacity_;
}

template <typename T, class A>
inline typename Buffer<T, A>::Pointer Buffer<T, A>::Data() QUOTIENT_NOEXCEPT {
  return data_;
}

template <typename T, class A>
inline typename Buffer<T, A>::ConstPointer Buffer<T, A>::Data() const
    QUOTIENT_NOEXCEPT {
  return data_;
}

template <typename T, class A>
inline typename Buffer<T, A>::Pointer Buffer<T, A>::begin() QUOTIENT_NOEXCEPT {
  return data_;
}

template <typename T, class A>
inline typename Buffer<T, A>::ConstPointer Buffer<T, A>::begin() const
    QUOTIENT_NOEXCEPT {
  return data_;
}

template <typename T, class A>
inline typename Buffer<T, A>::ConstPointer Buffer<T, A>::cbegin() const
    QUOTIENT_NOEXCEPT {
  return data_;
}

template <typename T, class A>
inline typename Buffer<T, A>::Pointer Buffer<T, A>::end() QUOTIENT_NOEXCEPT {
  return data_ + size_;
}

template <typename T, class A>
inline typename Buffer<T, A>::ConstPointer Buffer<T, A>::end() const
    QUOTIENT_NOEXCEPT {
  return data_ + size_;
}

template <typename T, class A>
inline typename Buffer<T, A>::ConstPointer Buffer<T, A>::cend() const
    QUOTIENT_NOEXCEPT {
  return data_ + size_;
}

template <typename T, class A>
inline typename Buffer<T, A>::Reference Buffer<T, A>::operator[](SizeType index)
    QUOTIENT_NOEXCEPT {
  return data_[index];
}

template <typename T, class A>
inline typename Buffer<T, A>::ConstReference Buffer<T, A>::operator[](
    SizeType index) const QUOTIENT_NOEXCEPT {
  return data_[index];
}

template <typename T, class A>
inline typename Buffer<T, A>::Reference Buffer<T, A>::Back() QUOTIENT_NOEXCEPT {
  return data_[size_ - 1];
}

template <typename T, class A>
inline typename Buffer<T, A>::ConstReference Buffer<T, A>::Back() const
    QUOTIENT_NOEXCEPT {
  return data_[size_ - 1];
}

template <typename T, class A>
inline void Buffer<T, A>::Clear() QUOTIENT_NOEXCEPT {
  DestructData();
  size_ = 0;
  capacity_ = 0;
  data_ = nullptr;
}

}  // namespace quotient

#endif  // ifndef QUOTIENT_BUFFER_IMPL_H_
