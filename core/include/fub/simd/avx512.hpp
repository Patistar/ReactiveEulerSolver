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

#ifndef AVX512_HPP
#define AVX512_HPP

#include "fub/simd/scalar.hpp"

#if __has_include("immintrin.h") && defined(__AVX512F__)
#define FUB_HAS_AVX512 1
#include "immintrin.h"

namespace fub {

/// @brief Simd wrapper for the SSE case (SISD case).
template <> class Simd<double, simd_abi::AVX512F> {
public:
  using value_type = double;
  using abi_type = simd_abi::AVX512F;
  using vector_type = __m512d;
  static constexpr int size() noexcept { return 8; }

  Simd() = default;

  // Implicit value conversion
  template <typename U,
            std::enable_if_t<std::is_convertible_v<U, value_type>>* = nullptr>
  Simd(U&& u) : m_value{_mm512_set1_pd(u)} {}

  // Implicit vector cast
  template <typename U,
            std::enable_if_t<std::is_convertible_v<U, vector_type>>* = nullptr>
  Simd(U&& u) : m_value{u} {}

  // Copy members

  void copy_from(const double* mem) noexcept { m_value = _mm512_load_pd(mem); }
  void copy_to(double* mem) const noexcept { _mm512_store_pd(mem, m_value); }

  friend Simd operator+(const Simd& x, const Simd& y) noexcept {
    return _mm512_add_pd(x.m_value, y.m_value);
  }

  friend Simd operator-(const Simd& x, const Simd& y) noexcept {
    return _mm512_min_pd(x.m_value, y.m_value);
  }

  friend Simd operator/(const Simd& x, const Simd& y) noexcept {
    return _mm512_div_pd(x.m_value, y.m_value);
  }

  friend Simd operator*(const Simd& x, const Simd& y) noexcept {
    return _mm512_mul_pd(x.m_value, y.m_value);
  }

  friend Simd sqrt(const Simd& x) noexcept {
    return _mm512_sqrt_pd(x.m_value);
  }

private:
  vector_type m_value;
};

} // namespace fub

#endif
#endif // !AVX512_HPP
