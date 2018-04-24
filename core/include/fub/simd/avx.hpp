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

#ifndef FUB_SIMD_AVX_HPP
#define FUB_SIMD_AVX_HPP

#include "fub/simd/scalar.hpp"

#if __has_include("immintrin.h") && defined(__AVX__) && !defined(FUB_NO_AVX)
#define FUB_SIMD_HAS_AVX
#include "immintrin.h"
#include <boost/align/is_aligned.hpp>

namespace fub {
namespace simd_abi {
struct avx {};
} // namespace simd_abi

template <typename T> using avx = simd<T, simd_abi::avx>;

template <> struct is_simd_abi<simd_abi::avx> : std::true_type {};
template <> struct is_simd<simd<double, simd_abi::avx>> : std::true_type {};
template <>
struct is_simd_mask<simd_mask<double, simd_abi::avx>> : std::true_type {};
template <typename U>
struct memory_alignment<simd<double, simd_abi::avx>, U>
    : std::integral_constant<int, 32> {};

template <> class simd_mask<double, simd_abi::avx> {
public:
  using value_type = bool;
  using abi_type = simd_abi::avx;
  using simd_type = simd<double, abi_type>;
  using native_type = __m256i;

  simd_mask() = default;

  explicit simd_mask(bool value) noexcept
      : m_value{value ? _mm256_set1_epi64x(0xFFFFFFFFFFFFFFFFul)
                      : _mm256_setzero_si256()} {}

  explicit simd_mask(const native_type& value) noexcept : m_value{value} {}
  explicit operator native_type() const noexcept { return m_value; }

  // Load and Store from Memory

  template <typename U>
  std::enable_if_t<std::is_convertible<U, value_type>::value, void>
  copy_from(const U* mem, ElementAlignmentTag) noexcept {
    m_value = _mm256_loadu_si256(mem);
  }

  template <typename U>
  std::enable_if_t<std::is_convertible<U, value_type>::value, void>
  copy_from(const U* mem, VectorAlignmentTag) noexcept {
    assert(boost::alignment::is_aligned(memory_alignment_v<simd_type, U>, mem));
    m_value = _mm256_load_si256(mem);
  }

  template <typename U>
  std::enable_if_t<std::is_convertible<U, value_type>::value, void>
  copy_to(U* mem, VectorAlignmentTag) const noexcept {
    assert(boost::alignment::is_aligned(memory_alignment_v<simd_type, U>, mem));
    _mm256_store_si256(mem, m_value);
  }

  template <typename U>
  std::enable_if_t<std::is_convertible<U, value_type>::value, void>
  copy_to(U* mem, ElementAlignmentTag) const noexcept {
    _mm256_storeu_si256(mem, m_value);
  }

  // Element Access

  // Unary Operators (element-wise)

  simd_mask operator!() const noexcept {
    return simd_mask{
        _mm256_xor_si256(m_value, _mm256_cmpeq_epi64(m_value, m_value))};
  }

  // Binary Operators (element-wise)

  friend simd_mask operator&&(const simd_mask& m1,
                              const simd_mask& m2) noexcept {
    return simd_mask{_mm256_and_si256(m1.m_value, m2.m_value)};
  }

  friend simd_mask operator||(const simd_mask& m1,
                              const simd_mask& m2) noexcept {
    return simd_mask{_mm256_or_si256(m1.m_value, m2.m_value)};
  }

  // Compare Operators (element-wise)

  friend simd_mask operator==(const simd_mask& m1,
                              const simd_mask& m2) noexcept {
    return simd_mask{_mm256_cmpeq_epi64(m1.m_value, m2.m_value)};
  }

  friend simd_mask operator!=(const simd_mask& m1,
                              const simd_mask& m2) noexcept {
    return !(m1 == m2);
  }

private:
  native_type m_value;
};

inline bool any_of(const simd_mask<double, simd_abi::avx>& x) noexcept {
  __m256i zero = _mm256_setzero_si256();
  __m256i truei = _mm256_cmpeq_epi8(zero, zero);
  __m256i v = static_cast<__m256i>(x);
  __m256i any = _mm256_cmpeq_epi8(v, truei);
  int mask = _mm256_movemask_epi8(any);
  return (mask != 0);
}

inline bool all_of(const simd_mask<double, simd_abi::avx>& x) noexcept {
  __m256i v = static_cast<__m256i>(x);
  __m256i vcmp = _mm256_cmpeq_epi8(v, _mm256_setzero_si256());
  int mask = _mm256_movemask_epi8(vcmp);
  return (mask == 0);
}

/// @brief simd wrapper for the avx case (SISD case).
template <> class simd<double, simd_abi::avx> {
public:
  using value_type = double;
  using abi_type = simd_abi::avx;
  using mask_type = simd_mask<value_type, abi_type>;
  using native_type = __m256d;
  static constexpr int size() noexcept { return 4; }

  simd() = default;

  // Implicit value broadcast
  template <
      typename U,
      std::enable_if_t<std::is_convertible<U, value_type>::value>* = nullptr>
  simd(U&& u) : m_value{_mm256_set1_pd(u)} {}

  // Implicit cast from native type
  simd(const native_type& value) noexcept : m_value{value} {}
  // Explicit cast into native type
  explicit operator const native_type&() const noexcept { return m_value; }
  explicit operator native_type&() noexcept { return m_value; }

  // Copy members

  template <typename U>
  std::enable_if_t<std::is_convertible<U, value_type>::value, void>
  copy_from(const U* mem, ElementAlignmentTag) noexcept {
    m_value = _mm256_loadu_pd(mem);
  }

  template <typename U>
  std::enable_if_t<std::is_convertible<U, value_type>::value, void>
  copy_to(U* mem, ElementAlignmentTag) const noexcept {
    _mm256_storeu_pd(mem, m_value);
  }

  template <typename U>
  std::enable_if_t<std::is_convertible<U, value_type>::value, void>
  copy_from(const U* mem, VectorAlignmentTag) noexcept {
    assert(boost::alignment::is_aligned(memory_alignment_v<simd, U>, mem));
    m_value = _mm256_load_pd(mem);
  }

  template <typename U>
  std::enable_if_t<std::is_convertible<U, value_type>::value, void>
  copy_to(U* mem, VectorAlignmentTag) const noexcept {
    assert(boost::alignment::is_aligned(memory_alignment_v<simd, U>, mem));
    _mm256_store_pd(mem, m_value);
  }

  friend simd operator+(const simd& x, const simd& y) noexcept {
    return _mm256_add_pd(x.m_value, y.m_value);
  }

  friend simd operator-(const simd& x, const simd& y) noexcept {
    return _mm256_sub_pd(x.m_value, y.m_value);
  }

  friend simd operator/(const simd& x, const simd& y) noexcept {
    return _mm256_div_pd(x.m_value, y.m_value);
  }

  friend simd operator*(const simd& x, const simd& y) noexcept {
    return _mm256_mul_pd(x.m_value, y.m_value);
  }

  simd operator-() const noexcept { return (*this) * (-1); }

  simd& operator+=(const simd& x) noexcept {
    m_value = _mm256_add_pd(m_value, x.m_value);
    return *this;
  }

  simd& operator-=(const simd& x) noexcept {
    m_value = _mm256_sub_pd(m_value, x.m_value);
    return *this;
  }

  simd& operator*=(const simd& x) noexcept {
    m_value = _mm256_mul_pd(m_value, x.m_value);
    return *this;
  }

  simd& operator/=(const simd& x) noexcept {
    m_value = _mm256_div_pd(m_value, x.m_value);
    return *this;
  }

  friend mask_type operator<(const simd& x, const simd& y) noexcept {
    return mask_type{_mm256_cmp_pd(x.m_value, y.m_value, _CMP_LT_OS)};
  }
  friend mask_type operator<=(const simd& x, const simd& y) noexcept {
    return mask_type{_mm256_cmp_pd(x.m_value, y.m_value, _CMP_LE_OS)};
  }
  friend mask_type operator>(const simd& x, const simd& y) noexcept {
    return mask_type{_mm256_cmp_pd(x.m_value, y.m_value, _CMP_GT_OS)};
  }
  friend mask_type operator>=(const simd& x, const simd& y) noexcept {
    return mask_type{_mm256_cmp_pd(x.m_value, y.m_value, _CMP_GE_OS)};
  }

  friend mask_type operator==(const simd& x, const simd& y) noexcept {
    return mask_type{_mm256_cmp_pd(x.m_value, y.m_value, _CMP_EQ_OS)};
  }
  friend mask_type operator!=(const simd& x, const simd& y) noexcept {
    return mask_type{_mm256_cmp_pd(x.m_value, y.m_value, _CMP_NEQ_OS)};
  }

private:
  native_type m_value;
};

template <>
class const_where_expression<simd_mask<double, simd_abi::avx>,
                             simd<double, simd_abi::avx>> {
public:
  using value_type = double;
  using abi_type = simd_abi::avx;
  using mask_type = simd_mask<double, simd_abi::avx>;
  using simd_type = simd<double, simd_abi::avx>;
  using native_type = typename simd_type::native_type;

  const_where_expression(const mask_type& m, const simd_type& s) noexcept
      : m_mask{m}, m_simd{const_cast<simd_type&>(s)} {}

  const_where_expression(const mask_type& m, simd_type& s) noexcept
      : m_mask{m}, m_simd{s} {}

#if __cplusplus >= 201703
  const_where_expression(const const_where_expression&) = delete;
  const_where_expression& operator=(const const_where_expression&) = delete;
#else
  const_where_expression(const const_where_expression&);
  const_where_expression& operator=(const const_where_expression&);
#endif

  // simd_type operator-() const&& { return m_mask ? -m_simd : m_simd; }

  template <typename U>
  std::enable_if_t<std::is_convertible<double, U>::value, void>
  copy_to(U* mem, ElementAlignmentTag) const&& noexcept {
    _mm256_maskstore_pd(mem, m_mask, m_simd);
  }

  template <typename U>
  std::enable_if_t<std::is_convertible<double, U>::value, void>
  copy_to(U* mem, VectorAlignmentTag) const&& noexcept {
    _mm256_maskstore_pd(mem, m_mask, m_simd);
  }

private:
  const mask_type& m_mask;
  simd_type& m_simd;

  friend class where_expression<mask_type, simd_type>;
};

template <>
class where_expression<simd_mask<double, simd_abi::avx>,
                       simd<double, simd_abi::avx>>
    : public const_where_expression<simd_mask<double, simd_abi::avx>,
                                    simd<double, simd_abi::avx>> {
public:
  using value_type = double;
  using abi_type = simd_abi::avx;
  using mask_type = simd_mask<double, simd_abi::avx>;
  using simd_type = simd<double, simd_abi::avx>;
  using native_type = simd_type::native_type;

private:
  using Base = const_where_expression<mask_type, simd_type>;
  using Base::m_mask;
  using Base::m_simd;

public:
  where_expression(const mask_type& m, simd_type& s) : Base(m, s) {}

  where_expression(const where_expression&);
  where_expression& operator=(const where_expression&);

  template <typename U>
  std::enable_if_t<std::is_constructible<simd_type, U>::value, void>
  operator=(U&& u) {
    simd_type source{u};
    // m_simd = (m_mask & source) | (!m_mask & m_simd);
    m_simd = _mm256_or_pd(
        _mm256_and_pd(_mm256_castsi256_pd(static_cast<__m256i>(m_mask)),
                      static_cast<native_type>(source)),
        _mm256_andnot_pd(_mm256_castsi256_pd(static_cast<__m256i>(m_mask)),
                         static_cast<native_type>(m_simd)));
  }

  template <typename U, typename Flags>
      std::enable_if_t<is_simd_flag_type_v<Flags> &&
                           std::is_convertible<U, double>::value,
                       void>
      copy_from(const U* mem, Flags) && noexcept {
    simd_type source{_mm256_maskload_pd(mem, m_mask)};
    *this = source;
  }
};

////////////////////////////////////////////////////////////////////////////////
// Math Function Overloads

inline simd<double, simd_abi::avx>
sqrt(const simd<double, simd_abi::avx>& x) noexcept {
  return _mm256_sqrt_pd(static_cast<__m256d>(x));
}

inline simd<double, simd_abi::avx>
abs(const simd<double, simd_abi::avx>& x) noexcept {
  __m256i minus1i = _mm256_set1_epi32(-1);
  __m256d minus1d = _mm256_castsi256_pd(_mm256_srli_epi32(minus1i, 1));
  __m256d result = _mm256_and_pd(minus1d, static_cast<__m256d>(x));
  return result;
}

////////////////////////////////////////////////////////////////////////////////
// Special Math Function Overloads

inline simd<double, simd_abi::avx>
min(const simd<double, simd_abi::avx>& x,
    const simd<double, simd_abi::avx>& y) noexcept {
  return _mm256_min_pd(static_cast<__m256d>(x), static_cast<__m256d>(y));
}

inline simd<double, simd_abi::avx>
max(const simd<double, simd_abi::avx>& x,
    const simd<double, simd_abi::avx>& y) noexcept {
  return _mm256_max_pd(static_cast<__m256d>(x), static_cast<__m256d>(y));
}

} // namespace fub

#endif
#endif // !FUB_SIMD_AVX_HPP
