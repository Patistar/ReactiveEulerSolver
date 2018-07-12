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

#ifndef FUB_CORE_ACCESSOR_SIMD_HPP
#define FUB_CORE_ACCESSOR_SIMD_HPP

#include "fub/simd.hpp"

#ifndef NDEBUG
#include <boost/alignment/is_aligned.hpp>
#include <assert>
#endif

namespace fub {
template <typename T, typename Abi, typename Alignment> class simd_proxy {
public:
  using element_type = T;
  using pointer = T*;
  using simd_type = simd<T, Abi>;

  simd_proxy(pointer p) noexcept : m_pointer{p} {}
  simd_proxy(const simd_proxy&) = delete;
  simd_proxy& operator=(const simd_proxy&) = delete;
  simd_proxy(simd_proxy&&) = default;
  simd_proxy& operator=(simd_proxy&&) = default;

  operator simd_type() const noexcept {
    return simd_type{m_pointer, Alignment};
  }

  simd_proxy& operator=(const simd_type& v) noexcept {
    v.copy_to(m_pointer, Alignment());
  }

private:
  pointer m_pointer;
};

template <typename T, typename Abi> struct accessor_simd_aligned {
  using pointer = T*;
  using reference = simd_proxy<T, Abi, vector_alignment_tag>;
  static constexpr reference access(pointer origin, std::ptrdiff_t offset) noexcept {
    pointer ptr = origin + offset;
    assert(boost::alignment::is_aligned(ptr, memory_alignment_v<simd<T, Abi>>));
    return reference(ptr);
  }
  template <typename S> using rebind = accessor_simd_aligned<S, Abi>;
};

template <typename T, typename Abi> struct accessor_simd_unaligned {
  using pointer = T*;
  using reference = simd_proxy<T, Abi, element_alignment_tag>;
  static constexpr reference access(pointer origin, std::ptrdiff_t offset) noexcept {
    return reference(origin + offset);
  }
  template <typename S> using rebind = accessor_simd_unaligned<S, Abi>;
};

} // namespace fub

#endif
