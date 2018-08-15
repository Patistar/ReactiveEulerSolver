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

#ifndef FUB_CORE_ADDRESS_NATIVE_HPP
#define FUB_CORE_ADDRESS_NATIVE_HPP

#include <cstddef>

namespace fub {
inline namespace v1 {

struct address_native {
  using difference_type = std::ptrdiff_t;

  void* address{nullptr};

  constexpr address_native() = default;

  template <typename T>
  constexpr explicit address_native(T* ptr)
      : address{static_cast<void*>(ptr)} {}

  template <typename T>
  constexpr std::ptrdiff_t distance_as(address_native a) const noexcept {
    return static_cast<T*>(a.address) - static_cast<T*>(address);
  }

  template <typename T>
  constexpr address_native& advance_as(std::ptrdiff_t n) noexcept {
    address = static_cast<void*>(static_cast<T*>(address) + n);
    return *this;
  }

  constexpr explicit operator void*() const noexcept { return address; }
  constexpr explicit operator const void*() const noexcept { return address; }
};

} // namespace v1
} // namespace fub

#endif
