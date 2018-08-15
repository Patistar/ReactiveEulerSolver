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

#ifndef FUB_CORE_SIMD_ADDRESS_HPP
#define FUB_CORE_SIMD_ADDRESS_HPP

#include "fub/simd.hpp"

#include <boost/align/is_aligned.hpp>

#include <cassert>

namespace fub {
inline namespace v1 {

template <int Alignment> struct aligned_address {
  using difference_type = std::ptrdiff_t;

  void* address{nullptr};

  constexpr aligned_address() = default;

  constexpr int alignment() const noexcept { return Alignment; }

  constexpr aligned_address(void* ptr) : address{ptr} {
    assert(boost::alignment::is_aligned(ptr, Alignment));
  }

  constexpr aligned_address(const void* ptr) : address{const_cast<void*>(ptr)} {
    assert(boost::alignment::is_aligned(ptr, Alignment));
  }

  template <typename T>
  constexpr std::ptrdiff_t distance_as(aligned_address a) const noexcept {
    using S = remove_cvref_t<T>;
    constexpr int width = alignof(S) < Alignment ? Alignment / alignof(S) : 1;
    return (static_cast<S*>(a.address) - static_cast<S*>(address)) / width;
  }

  template <typename T>
  constexpr aligned_address& advance_as(std::ptrdiff_t n) noexcept {
    using S = remove_cvref_t<T>;
    constexpr int width = alignof(S) < Alignment ? Alignment / alignof(S)  : 1;
    address = static_cast<void*>(static_cast<S*>(address) + n * width);
    return *this;
  }

  constexpr bool is_null() const noexcept { return address != nullptr; }

  constexpr operator void*() const noexcept { return address; }
  constexpr operator const void*() const noexcept { return address; }
};

} // namespace v1
} // namespace fub

#endif