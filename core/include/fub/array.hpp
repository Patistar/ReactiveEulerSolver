// Copyright (c) 2018 Maikel Nadolski
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

/// @file This module provides a constexpr-enabled array type. The Standard
/// Library supports constexpr for arrays only from C++17 on upwards. In case
/// of an incompatible standard library implementation we provide one on our
/// own.

#ifndef FUB_CORE_ARRAY_HPP
#define FUB_CORE_ARRAY_HPP

#include "fub/core/config.hpp"

#ifdef FUB_CORE_USE_STD_ARRAY
#include <array>
namespace fub {
using std::array;
}
#else
#include <cstdint>
#include <range/v3/algorithm/equal.hpp>

namespace fub {
template <typename T, std::size_t N> struct array {
  using value_type = T;
  using pointer = T*;
  using const_pointer = const T*;
  using reference = T&;
  using const_reference = const T&;
  using iterator = pointer;
  using const_iterator = const_pointer;
  using difference_type = std::ptrdiff_t;
  using size_type = std::ptrdiff_t;

  T m_array[N];

  constexpr pointer data() noexcept { return &m_array[0]; }

  constexpr const_pointer data() const noexcept { return &m_array[0]; }

  constexpr const_pointer cdata() const noexcept { return &m_array[0]; }

  constexpr size_type size() const noexcept { return N; }

  constexpr reference operator[](difference_type i) noexcept {
    return m_array[i];
  }
  constexpr const_reference operator[](difference_type i) const noexcept {
    return m_array[i];
  }

  constexpr iterator begin() noexcept { return data(); }
  constexpr const_iterator begin() const noexcept { return data(); }
  constexpr const_iterator cbegin() const noexcept { return data(); }

  constexpr iterator end() noexcept { return data() + N; }
  constexpr const_iterator end() const noexcept { return data() + N; }
  constexpr const_iterator cend() const noexcept { return data() + N; }
};

template <typename T, std::size_t N>
constexpr bool operator==(const array<T, N>& x, const array<T, N>& y) noexcept {
  return ranges::equal(x, y);
}

template <typename T, std::size_t N>
constexpr bool operator!=(const array<T, N>& x, const array<T, N>& y) noexcept {
  return !(x == y);
}
} // namespace fub
#endif
#endif
