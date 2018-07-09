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
#include <fmt/ostream.h>
#include <ostream>

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

namespace detail {
template <typename T, typename Abi, int... Is>
simd<T, Abi> power_impl(std::integer_sequence<int, Is...>,
                        simd<T, Abi> x) noexcept {
  simd<T, Abi> one{1.0};
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-value"
#endif
  return (one * ... * (Is, x));
#ifdef __clang__
#pragma clang diagnostic pop
#endif
}

template <int N, typename T, typename Abi>
simd<T, Abi> power(simd<T, Abi> x) noexcept {
  return power_impl(std::make_integer_sequence<int, N>(), x);
}

template <typename T, typename Abi, typename... Ts, int... Is>
simd<T, Abi> polynomial_impl(std::integer_sequence<int, Is...>, simd<T, Abi> x,
                             Ts... cs) noexcept {
  return (x + ... + (power<Is + 1>(x) * cs));
}

template <typename T, typename Abi, typename... Ts>
simd<T, Abi> polynomial(simd<T, Abi> x, Ts... cs) noexcept {
  return polynomial_impl(std::make_integer_sequence<int, sizeof...(Ts)>(), x,
                         cs...);
}
} // namespace detail

template <typename T, typename Abi,
          typename = std::enable_if_t<std::is_floating_point<T>::value>>
simd<T, Abi> exp(simd<T, Abi> x) noexcept {
  // 1. / log(2)
  static constexpr T log2e = 1.44269504088896340736;
  // Taylor coefficients, 1/n!
  // Not using minimax approximation because we prioritize precision close
  // to x = 0
  static constexpr T p2 = 1. / 2.;
  static constexpr T p3 = 1. / 6.;
  static constexpr T p4 = 1. / 24.;
  static constexpr T p5 = 1. / 120.;
  static constexpr T p6 = 1. / 720.;
  static constexpr T p7 = 1. / 5040.;
  static constexpr T p8 = 1. / 40320.;
  static constexpr T p9 = 1. / 362880.;
  static constexpr T p10 = 1. / 3628800.;
  static constexpr T p11 = 1. / 39916800.;
  static constexpr T p12 = 1. / 479001600.;
  static constexpr T p13 = 1. / 6227020800.;
  static const T ln2d_hi = 0.693145751953125;
  static const T ln2d_lo = 1.42860682030941723212E-6;
  const simd<T, Abi> r = round(x * log2e);
  const simd<T, Abi> y = x - r * (ln2d_lo + ln2d_hi);
  return detail::polynomial(y, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12,
                            p13);
}

template <typename T, typename Abi,
          typename = std::enable_if_t<std::is_floating_point<T>::value>>
simd<T, Abi> log(simd<T, Abi> x) noexcept {
  const simd<T, Abi> z = (x - 1) / (x + 1);
  const simd<T, Abi> z3 = detail::power<3>(z) / 3.0;
  const simd<T, Abi> z5 = detail::power<5>(z) / 5.0;
  const simd<T, Abi> z7 = detail::power<7>(z) / 7.0;
  where(x > 0, x) = 2.0 * (z + z3 + z5 + z7);
  const simd<T, Abi> not_a_number = std::numeric_limits<T>::quiet_NaN();
  where(x <= 0, x) = not_a_number;
  return x;
}

template <typename T, typename Abi, typename = std::enable_if_t<std::is_floating_point<T>::value>>
simd<T, Abi> log10(simd<T, Abi> x) noexcept {
  const simd<T, Abi> ln10 = 2.30258509299;
  return ln10 * log(x);
}

template <typename T, typename Abi,
          typename = std::enable_if_t<std::is_floating_point<T>::value>>
simd<T, Abi> pow(simd<T, Abi> a, simd<T, Abi> b) noexcept {
  return exp(a * log(b));
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
