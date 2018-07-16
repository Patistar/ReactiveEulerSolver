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

#ifndef FUB_CORE_ACCESSOR_NATIVE_HPP
#define FUB_CORE_ACCESSOR_NATIVE_HPP

#include <cstddef>

namespace fub {
template <typename T = void> struct accessor_native {
  using pointer = T*;
  using reference = T&;
  static constexpr reference access(pointer ptr, std::ptrdiff_t /* size */,
                                    std::ptrdiff_t n) noexcept {
    return *(ptr + n);
  }
  template <typename S> using rebind = accessor_native<S>;
};

template <> struct accessor_native<void> {
  using pointer = void*;
  template <typename S> using rebind = accessor_native<S>;
};
} // namespace fub

#endif
