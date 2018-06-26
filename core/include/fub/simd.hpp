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

#ifndef FUB_CORE_SIMD_HPP
#define FUB_CORE_SIMD_HPP

#include "fub/simd/avx.hpp"
#include "fub/simd/avx512.hpp"
#include "fub/simd/native.hpp"
#include "fub/simd/scalar.hpp"
#include "fub/simd/sse.hpp"

#include "fub/utility.hpp"

#include <array>
#include <cassert>
#include <ostream>
#include <fmt/ostream.h>

namespace fub {

template <typename T, typename Abi>
simd<T, Abi> clamp(const simd<T, Abi>& v, const nodeduce_t<simd<T, Abi>>& lo,
                   const nodeduce_t<simd<T, Abi>>& hi) {
  assert(all_of(lo < hi));
  simd<T, Abi> x{v};
  where(x < lo, x) = lo;
  where(hi < x, x) = hi;
  return x;
}

template <typename A, typename Abi, typename I>
simd_mask<A, Abi> almost_equal(const simd<A, Abi>& x,
                               const nodeduce_t<simd<A, Abi>>& y,
                               const I& ulp) noexcept {
  constexpr A eps = std::numeric_limits<A>::epsilon();
  constexpr A min = std::numeric_limits<A>::min();
  const simd<A, Abi> diff = fub::abs(x - y);
  const simd<A, Abi> sum = fub::abs(x + y);
  return diff <= (eps * sum * ulp) || diff < min;
}

template <typename T, typename Abi>
std::ostream& operator<<(std::ostream& out, const simd<T, Abi>& v) {
  constexpr int width = simd_size_v<simd<T, Abi>>;
  constexpr int alignment = memory_alignment_v<simd<T, Abi>, T>;
  alignas(alignment) std::array<T, width> array;
  v.copy_to(array.data(), vector_alignment);
  out << "(" << array[0];
  for (std::size_t i = 1; i < width; ++i) {
    out << ", " << array[i];
  }
  out << ")";
  return out;
}

} // namespace fub

#endif
