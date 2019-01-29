/*
 * Copyright (c) 2019 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_BUFFER_COMPLEX_IMPL_H_
#define QUOTIENT_BUFFER_COMPLEX_IMPL_H_

#include <algorithm>

#include "quotient/buffer.hpp"

namespace quotient {

template <typename Real>
inline Buffer<Complex<Real>>::Buffer() noexcept : size_(0),
                                                  capacity_(0),
                                                  data_(nullptr) {}

template <typename Real>
inline void Buffer<Complex<Real>>::DestructData() {
  Real* real_data = reinterpret_cast<Real*>(data_);

  if (!is_trivially_destructible_real) {
    const SizeType real_size = 2 * size_;
    for (Real* iter = real_data; iter != real_data + real_size; ++iter) {
      iter->~Real();
    }
  }

  const SizeType real_capacity = 2 * capacity_;
  AllocatorTraits::deallocate(allocator_, real_data, real_capacity);
}

template <typename Real>
inline Buffer<Complex<Real>>::~Buffer() {
  DestructData();
}

template <typename Real>
inline void Buffer<Complex<Real>>::ConstructRange(SizeType start,
                                                  SizeType end) {
  if (!is_trivially_constructible_real) {
    Real* real_data = reinterpret_cast<Real*>(data_);
    const SizeType real_start = 2 * start;
    const SizeType real_end = 2 * end;
    for (Real* iter = real_data + real_start; iter != real_data + real_end;
         ++iter) {
      AllocatorTraits::construct(allocator_, iter);
    }
  }
}

template <typename Real>
inline Buffer<Complex<Real>>::Buffer(SizeType num_elements)
    : size_(num_elements), capacity_(num_elements) {
  const SizeType real_capacity = 2 * capacity_;
  data_ = reinterpret_cast<Complex<Real>*>(
      AllocatorTraits::allocate(allocator_, real_capacity));
  ConstructRange(0, size_);
}

template <typename Real>
inline void Buffer<Complex<Real>>::FillConstructRange(SizeType start,
                                                      SizeType end,
                                                      ConstReference value) {
  if (is_trivially_copy_constructible_real) {
    std::fill(data_ + start, data_ + end, value);
  } else {
    Real* real_data = reinterpret_cast<Real*>(data_);
    for (SizeType offset = start; offset < end; ++offset) {
      AllocatorTraits::construct(allocator_, real_data + 2 * offset,
                                 value.real());
      AllocatorTraits::construct(allocator_, real_data + 2 * offset + 1,
                                 value.imag());
    }
  }
}

template <typename Real>
inline void Buffer<Complex<Real>>::CopyConstructRange(SizeType start,
                                                      ConstIterator begin,
                                                      ConstIterator end) {
  if (is_trivially_copy_constructible_real) {
    std::copy(begin, end, data_ + start);
  } else {
    Real* real_data = reinterpret_cast<Real*>(data_);
    SizeType offset = start;
    for (ConstIterator iter = begin; iter != end; ++iter, ++offset) {
      AllocatorTraits::construct(allocator_, real_data + 2 * offset,
                                 iter->real());
      AllocatorTraits::construct(allocator_, real_data + 2 * offset + 1,
                                 iter->imag());
    }
  }
}

template <typename Real>
inline Buffer<Complex<Real>>::Buffer(SizeType num_elements,
                                     ConstReference value)
    : size_(num_elements), capacity_(num_elements) {
  const SizeType real_capacity = 2 * capacity_;
  data_ = reinterpret_cast<Complex<Real>*>(
      AllocatorTraits::allocate(allocator_, real_capacity));
  FillConstructRange(0, size_, value);
}

template <typename Real>
Buffer<Complex<Real>>::Buffer(ConstIterator begin, ConstIterator end) {
  size_ = std::distance(begin, end);
  capacity_ = size_;

  const SizeType real_capacity = 2 * capacity_;
  data_ = reinterpret_cast<Complex<Real>*>(
      AllocatorTraits::allocate(allocator_, real_capacity));

  CopyConstructRange(0, begin, end);
}

template <typename Real>
Buffer<Complex<Real>>::Buffer(std::initializer_list<Complex<Real>> list)
    : Buffer(list.begin(), list.end()) {}

template <typename Real>
inline Buffer<Complex<Real>>::Buffer(const Buffer<Complex<Real>>& buffer)
    : Buffer(buffer.begin(), buffer.end()) {}

template <typename Real>
inline Buffer<Complex<Real>>::Buffer(const std::vector<Complex<Real>>& vec)
    : Buffer(vec.begin(), vec.end()) {}

template <typename Real>
inline Buffer<Complex<Real>>::Buffer(Buffer<Complex<Real>>&& buffer) noexcept
    : size_(buffer.size_),
      capacity_(buffer.capacity_),
      data_(buffer.data_) {
  buffer.size_ = 0;
  buffer.capacity_ = 0;
  buffer.data_ = nullptr;
}

template <typename Real>
Buffer<Complex<Real>>& Buffer<Complex<Real>>::operator=(
    const Buffer<Complex<Real>>& buffer) {
  if (this != &buffer) {
    const SizeType num_elements = buffer.Size();
    if (num_elements > capacity_) {
      DestructData();
      const SizeType real_capacity = 2 * num_elements;
      data_ = reinterpret_cast<Complex<Real>*>(
          AllocatorTraits::allocate(allocator_, real_capacity));
      size_ = num_elements;
      capacity_ = num_elements;
      CopyConstructRange(0, buffer.begin(), buffer.end());
    } else if (num_elements > size_) {
      std::copy(buffer.begin(), buffer.begin() + size_, data_);
      CopyConstructRange(size_, buffer.begin() + size_, buffer.end());
      size_ = num_elements;
    } else {
      std::copy(buffer.begin(), buffer.begin() + num_elements, data_);
      size_ = num_elements;
    }
  }
  return *this;
}

template <typename Real>
Buffer<Complex<Real>>& Buffer<Complex<Real>>::operator=(
    Buffer<Complex<Real>>&& buffer) noexcept {
  size_ = buffer.size_;
  capacity_ = buffer.capacity_;
  data_ = buffer.data_;

  buffer.size_ = 0;
  buffer.capacity_ = 0;
  buffer.data_ = nullptr;
}

template <typename Real>
Buffer<Complex<Real>>& Buffer<Complex<Real>>::operator=(
    const std::vector<Complex<Real>>& vec) {
  const SizeType num_elements = vec.size();
  if (num_elements > capacity_) {
    DestructData();
    const SizeType real_capacity = 2 * num_elements;
    data_ = reinterpret_cast<Complex<Real>*>(
        AllocatorTraits::allocate(allocator_, real_capacity));
    size_ = num_elements;
    capacity_ = num_elements;
    CopyConstructRange(0, vec.begin(), vec.end());
  } else if (num_elements > size_) {
    std::copy(vec.begin(), vec.begin() + size_, data_);
    CopyConstructRange(size_, vec.begin() + size_, vec.end());
    size_ = num_elements;
  } else {
    std::copy(vec.begin(), vec.begin() + num_elements, data_);
    size_ = num_elements;
  }

  return *this;
}

template <typename Real>
bool Buffer<Complex<Real>>::operator==(
    const Buffer<Complex<Real>>& buffer) const {
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

template <typename Real>
bool Buffer<Complex<Real>>::operator!=(
    const Buffer<Complex<Real>>& buffer) const {
  return !this->operator==(buffer);
}

template <typename Real>
void Buffer<Complex<Real>>::Resize(SizeType num_elements) {
  if (num_elements > capacity_) {
    DestructData();
    const SizeType real_capacity = 2 * num_elements;
    data_ = reinterpret_cast<Complex<Real>*>(
        AllocatorTraits::allocate(allocator_, real_capacity));
    size_ = num_elements;
    capacity_ = num_elements;
    ConstructRange(0, size_);
  } else if (num_elements > size_) {
    ConstructRange(size_, num_elements);
    size_ = num_elements;
  } else {
    size_ = num_elements;
  }
}

template <typename Real>
void Buffer<Complex<Real>>::Resize(SizeType num_elements,
                                   ConstReference value) {
  if (num_elements > capacity_) {
    DestructData();
    const SizeType real_capacity = 2 * num_elements;
    data_ = reinterpret_cast<Complex<Real>*>(
        AllocatorTraits::allocate(allocator_, real_capacity));
    size_ = num_elements;
    capacity_ = num_elements;
    FillConstructRange(0, size_, value);
  } else if (num_elements > size_) {
    std::fill(data_, data_ + size_, value);
    FillConstructRange(size_, num_elements, value);
    size_ = num_elements;
  } else {
    std::fill(data_, data_ + num_elements, value);
    size_ = num_elements;
  }
}

template <typename Real>
inline typename Buffer<Complex<Real>>::SizeType Buffer<Complex<Real>>::Size()
    const noexcept {
  return size_;
}

template <typename Real>
inline bool Buffer<Complex<Real>>::Empty() const noexcept {
  return size_ == 0;
}

template <typename Real>
inline typename Buffer<Complex<Real>>::SizeType
Buffer<Complex<Real>>::Capacity() const noexcept {
  return capacity_;
}

template <typename Real>
inline typename Buffer<Complex<Real>>::Pointer
Buffer<Complex<Real>>::Data() noexcept {
  return data_;
}

template <typename Real>
inline typename Buffer<Complex<Real>>::ConstPointer
Buffer<Complex<Real>>::Data() const noexcept {
  return data_;
}

template <typename Real>
inline typename Buffer<Complex<Real>>::Pointer
Buffer<Complex<Real>>::begin() noexcept {
  return data_;
}

template <typename Real>
inline typename Buffer<Complex<Real>>::ConstPointer
Buffer<Complex<Real>>::begin() const noexcept {
  return data_;
}

template <typename Real>
inline typename Buffer<Complex<Real>>::ConstPointer
Buffer<Complex<Real>>::cbegin() const noexcept {
  return data_;
}

template <typename Real>
inline typename Buffer<Complex<Real>>::Pointer
Buffer<Complex<Real>>::end() noexcept {
  return data_ + size_;
}

template <typename Real>
inline typename Buffer<Complex<Real>>::ConstPointer Buffer<Complex<Real>>::end()
    const noexcept {
  return data_ + size_;
}

template <typename Real>
inline typename Buffer<Complex<Real>>::ConstPointer
Buffer<Complex<Real>>::cend() const noexcept {
  return data_ + size_;
}

template <typename Real>
inline
    typename Buffer<Complex<Real>>::Reference Buffer<Complex<Real>>::operator[](
        SizeType index) noexcept {
  return data_[index];
}

template <typename Real>
inline typename Buffer<Complex<Real>>::ConstReference Buffer<Complex<Real>>::
operator[](SizeType index) const noexcept {
  return data_[index];
}

template <typename Real>
inline typename Buffer<Complex<Real>>::Reference
Buffer<Complex<Real>>::Back() noexcept {
  return data_[size_ - 1];
}

template <typename Real>
inline typename Buffer<Complex<Real>>::ConstReference
Buffer<Complex<Real>>::Back() const noexcept {
  return data_[size_ - 1];
}

template <typename Real>
inline void Buffer<Complex<Real>>::Clear() noexcept {
  DestructData();

  size_ = 0;
  capacity_ = 0;
  data_ = nullptr;
}

}  // namespace quotient

#endif  // ifndef QUOTIENT_BUFFER_COMPLEX_IMPL_H_
