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

#ifndef FUB_SIMD_SSE_HPP
#define FUB_SIMD_SSE_HPP

#include "fub/simd/scalar.hpp"

#if __has_include("emmintrin.h") && defined(__SSE__) && !defined(FUB_NO_SSE)
#define FUB_SIMD_HAS_SSE
#include "emmintrin.h"
#include "immintrin.h"
#include <boost/align/is_aligned.hpp>

namespace fub {
namespace simd_abi {
struct sse {};
} // namespace simd_abi

template <typename T> using sse = simd<T, simd_abi::sse>;

template <> struct is_simd_abi<simd_abi::sse> : std::true_type {};
template <> struct is_simd<simd<double, simd_abi::sse>> : std::true_type {};
template <>
struct is_simd_mask<simd_mask<double, simd_abi::sse>> : std::true_type {};
template <typename U>
struct memory_alignment<simd<double, simd_abi::sse>, U>
    : std::integral_constant<int, 16> {};

template <> class simd_mask<double, simd_abi::sse> {
public:
  using value_type = bool;
  using abi_type = simd_abi::sse;
  using simd_type = simd<double, abi_type>;
  using native_type = __m128d;

  simd_mask() = default;

  explicit simd_mask(bool value) noexcept
      : m_value{value ? _mm_set1_pd(0xFFFFFFFFFFFFFFFFul) : _mm_set1_pd(0ul)} {}

  explicit simd_mask(const native_type& value) noexcept : m_value{value} {}
  explicit operator native_type() const noexcept { return m_value; }

  // Load and Store from Memory

  template <typename U>
  std::enable_if_t<std::is_convertible<U, value_type>::value, void>
  copy_from(const U* mem, element_alignment_tag) noexcept {
    m_value = _mm_loadu_pd(mem);
  }

  template <typename U>
  std::enable_if_t<std::is_convertible<U, value_type>::value, void>
  copy_from(const U* mem, vector_alignment_tag) noexcept {
    assert(boost::alignment::is_aligned(memory_alignment_v<simd_type, U>, mem));
    m_value = _mm_load_pd(mem);
  }

  template <typename U>
  std::enable_if_t<std::is_convertible<U, value_type>::value, void>
  copy_to(U* mem, vector_alignment_tag) const noexcept {
    assert(boost::alignment::is_aligned(memory_alignment_v<simd_type, U>, mem));
    _mm_store_pd(mem, m_value);
  }

  template <typename U>
  std::enable_if_t<std::is_convertible<U, value_type>::value, void>
  copy_to(U* mem, element_alignment_tag) const noexcept {
    _mm_storeu_pd(mem, m_value);
  }

  // Element Access

  // Unary Operators (element-wise)

  simd_mask operator!() const noexcept {
    __m128i valuei = _mm_castpd_si128(m_value);
    __m128i truei = _mm_cmpeq_epi32(valuei, valuei);
    __m128i resulti = _mm_xor_si128(valuei, truei);
    __m128d resultd = _mm_castsi128_pd(resulti);
    return simd_mask{resultd};
  }

  // Binary Operators (element-wise)

  friend simd_mask operator&&(const simd_mask& m1,
                              const simd_mask& m2) noexcept {
    return simd_mask{_mm_and_pd(m1.m_value, m2.m_value)};
  }

  friend simd_mask operator||(const simd_mask& m1,
                              const simd_mask& m2) noexcept {
    return simd_mask{_mm_or_pd(m1.m_value, m2.m_value)};
  }

  // Compare Operators (element-wise)

  friend simd_mask operator==(const simd_mask& m1,
                              const simd_mask& m2) noexcept {
    return simd_mask{_mm_cmpeq_pd(m1.m_value, m2.m_value)};
  }

  friend simd_mask operator!=(const simd_mask& m1,
                              const simd_mask& m2) noexcept {
    return simd_mask{_mm_cmpneq_pd(m1.m_value, m2.m_value)};
  }

private:
  native_type m_value{};
};

inline bool any_of(const simd_mask<double, simd_abi::sse>& x) noexcept {
  __m128d zero = _mm_setzero_pd();
  __m128d trued = _mm_cmpeq_pd(zero, zero);
  __m128i truei = _mm_castpd_si128(trued);
  __m128i v = _mm_castpd_si128(static_cast<__m128d>(x));
  __m128i any = _mm_cmpeq_epi8(v, truei);
  int mask = _mm_movemask_epi8(any);
  return (mask != 0);
}

inline bool all_of(const simd_mask<double, simd_abi::sse>& x) noexcept {
  __m128i v = _mm_castpd_si128(static_cast<__m128d>(x));
  __m128i vcmp = _mm_cmpeq_epi8(v, _mm_setzero_si128());
  int mask = _mm_movemask_epi8(vcmp);
  return (mask == 0);
}

/// @brief simd wrapper for the sse case (SISD case).
template <> class simd<double, simd_abi::sse> {
public:
  using value_type = double;
  using abi_type = simd_abi::sse;
  using mask_type = simd_mask<value_type, abi_type>;
  using native_type = __m128d;
  static constexpr int size() noexcept { return 2; }

  simd() = default;

  // Implicit value broadcast
  template <
      typename U,
      std::enable_if_t<std::is_convertible<U, value_type>::value>* = nullptr>
  simd(U&& u) : m_value{_mm_set1_pd(u)} {}

  // Implicit cast from native type
  simd(const native_type& value) noexcept : m_value{value} {}
  // Explicit cast into native type
  explicit operator const native_type&() const noexcept { return m_value; }
  explicit operator native_type&() noexcept { return m_value; }

  // Copy members

  template <typename U>
  std::enable_if_t<std::is_convertible<U, value_type>::value, void>
  copy_from(const U* mem, element_alignment_tag) noexcept {
    m_value = _mm_loadu_pd(mem);
  }

  template <typename U>
  std::enable_if_t<std::is_convertible<U, value_type>::value, void>
  copy_to(U* mem, element_alignment_tag) const noexcept {
    _mm_storeu_pd(mem, m_value);
  }

  template <typename U>
  std::enable_if_t<std::is_convertible<U, value_type>::value, void>
  copy_from(const U* mem, vector_alignment_tag) noexcept {
    assert(boost::alignment::is_aligned(memory_alignment_v<simd, U>, mem));
    m_value = _mm_load_pd(mem);
  }

  template <typename U>
  std::enable_if_t<std::is_convertible<U, value_type>::value, void>
  copy_to(U* mem, vector_alignment_tag) const noexcept {
    assert(boost::alignment::is_aligned(memory_alignment_v<simd, U>, mem));
    _mm_store_pd(mem, m_value);
  }

  friend simd operator+(const simd& x, const simd& y) noexcept {
    return _mm_add_pd(x.m_value, y.m_value);
  }

  friend simd operator-(const simd& x, const simd& y) noexcept {
    return _mm_sub_pd(x.m_value, y.m_value);
  }

  friend simd operator/(const simd& x, const simd& y) noexcept {
    return _mm_div_pd(x.m_value, y.m_value);
  }

  friend simd operator*(const simd& x, const simd& y) noexcept {
    return _mm_mul_pd(x.m_value, y.m_value);
  }

  simd operator-() const noexcept { return (*this) * (-1); }

  simd& operator+=(const simd& x) noexcept {
    m_value = _mm_add_pd(m_value, x.m_value);
    return *this;
  }

  simd& operator-=(const simd& x) noexcept {
    m_value = _mm_sub_pd(m_value, x.m_value);
    return *this;
  }

  simd& operator*=(const simd& x) noexcept {
    m_value = _mm_mul_pd(m_value, x.m_value);
    return *this;
  }

  simd& operator/=(const simd& x) noexcept {
    m_value = _mm_div_pd(m_value, x.m_value);
    return *this;
  }

  friend mask_type operator<(const simd& x, const simd& y) noexcept {
    return mask_type{_mm_cmplt_pd(x.m_value, y.m_value)};
  }
  friend mask_type operator<=(const simd& x, const simd& y) noexcept {
    return mask_type{_mm_cmple_pd(x.m_value, y.m_value)};
  }
  friend mask_type operator>(const simd& x, const simd& y) noexcept {
    return mask_type{_mm_cmpgt_pd(x.m_value, y.m_value)};
  }
  friend mask_type operator>=(const simd& x, const simd& y) noexcept {
    return mask_type{_mm_cmpge_pd(x.m_value, y.m_value)};
  }

  friend mask_type operator==(const simd& x, const simd& y) noexcept {
    return mask_type{_mm_cmpeq_pd(x.m_value, y.m_value)};
  }
  friend mask_type operator!=(const simd& x, const simd& y) noexcept {
    return mask_type{_mm_cmpneq_pd(x.m_value, y.m_value)};
  }

private:
  native_type m_value{};
};

template <>
class const_where_expression<simd_mask<double, simd_abi::sse>,
                             simd<double, simd_abi::sse>> {
public:
  using value_type = double;
  using abi_type = simd_abi::sse;
  using mask_type = simd_mask<double, simd_abi::sse>;
  using simd_type = simd<double, simd_abi::sse>;
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
  copy_to(U* mem, element_alignment_tag) const&& noexcept {
    native_type value = _mm_and_pd(static_cast<native_type>(m_mask),
                                   static_cast<native_type>(m_simd));
    _mm_storeu_pd(mem, value);
  }

  template <typename U>
  std::enable_if_t<std::is_convertible<double, U>::value, void>
  copy_to(U* mem, vector_alignment_tag) const&& noexcept {
    native_type value = _mm_and_pd(static_cast<native_type>(m_mask),
                                   static_cast<native_type>(m_simd));
    assert(boost::alignment::is_aligned(memory_alignment_v<simd_type, U>, mem));
    _mm_store_pd(mem, value);
  }

private:
  const mask_type& m_mask;
  simd_type& m_simd;

  friend class where_expression<mask_type, simd_type>;
};

template <>
class where_expression<simd_mask<double, simd_abi::sse>,
                       simd<double, simd_abi::sse>>
    : public const_where_expression<simd_mask<double, simd_abi::sse>,
                                    simd<double, simd_abi::sse>> {
public:
  using value_type = double;
  using abi_type = simd_abi::sse;
  using mask_type = simd_mask<double, simd_abi::sse>;
  using simd_type = simd<double, simd_abi::sse>;
  using native_type = simd_type::native_type;

private:
  using Base = const_where_expression<mask_type, simd_type>;
  using Base::m_mask;
  using Base::m_simd;

public:
  where_expression(const mask_type& m, simd_type& s) : Base(m, s) {}

  where_expression(const where_expression&);
  where_expression& operator=(const where_expression&);

  template <typename U, typename Flags>
  std::enable_if_t<
      is_simd_flag_type_v<Flags> && std::is_convertible<U, double>::value, void>
  copy_from(const U* mem, Flags flags) const&& noexcept {
    simd_type source;
    source.copy_from(mem, flags);
    // m_simd = (m_mask & source) | (!m_mask & m_simd);
    m_simd = _mm_or_pd(_mm_and_pd(static_cast<native_type>(m_mask),
                                  static_cast<native_type>(source)),
                       _mm_andnot_pd(static_cast<native_type>(m_mask),
                                     static_cast<native_type>(m_simd)));
  }

  template <typename U>
  std::enable_if_t<std::is_constructible<simd_type, U>::value, void>
  operator=(U&& u) {
    simd_type source = u;
    m_simd = _mm_or_pd(_mm_and_pd(static_cast<native_type>(m_mask),
                                  static_cast<native_type>(source)),
                       _mm_andnot_pd(static_cast<native_type>(m_mask),
                                     static_cast<native_type>(m_simd)));
  }
};

////////////////////////////////////////////////////////////////////////////////
// Math Function Overloads

inline simd<double, simd_abi::sse>
sqrt(const simd<double, simd_abi::sse>& x) noexcept {
  return _mm_sqrt_pd(static_cast<__m128d>(x));
}

namespace sse_detail {
// Generate a constant vector of 4 integers stored in memory.
// Can be converted to any integer vector type
template <int i0, int i1, int i2, int i3> static inline __m128i constant4i() {
  static const union {
    int i[4];
    __m128i xmm;
  } u = {{i0, i1, i2, i3}};
  return u.xmm;
}
} // namespace sse_detail

inline simd<double, simd_abi::sse>
round(const simd<double, simd_abi::sse>& x) noexcept {
#ifdef __SSE4_2__
  return _mm_round_pd(static_cast<__m128d>(x),
                      _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
#else
  // Note: assume MXCSR control register is set to rounding
  // (don't use conversion to int, it will limit the value to +/- 2^31)
  __m128d signmask = _mm_castsi128_pd(
      sse_detail::constant4i<0, (int)0x80000000, 0, (int)0x80000000>()); // -0.0
  __m128d magic = _mm_castsi128_pd(
      sse_detail::constant4i<0, 0x43300000, 0,
                             0x43300000>()); // magic number = 2^52
  __m128d sign = _mm_and_pd(static_cast<__m128d>(x), signmask); // signbit of x
  sse<double> signedmagic = _mm_or_pd(magic, sign); // magic number with sign of x
  return x + signedmagic - signedmagic;         // round by adding magic number
#endif
}

inline simd<double, simd_abi::sse>
abs(const simd<double, simd_abi::sse>& x) noexcept {
  __m128i minus1i = _mm_set1_epi32(-1);
  __m128d minus1d = _mm_castsi128_pd(_mm_srli_epi32(minus1i, 1));
  __m128d result = _mm_and_pd(minus1d, static_cast<__m128d>(x));
  return result;
}

////////////////////////////////////////////////////////////////////////////////
// Special Math Function Overloads

inline simd<double, simd_abi::sse>
min(const simd<double, simd_abi::sse>& x,
    const simd<double, simd_abi::sse>& y) noexcept {
  return _mm_min_pd(static_cast<__m128d>(x), static_cast<__m128d>(y));
}

inline simd<double, simd_abi::sse>
max(const simd<double, simd_abi::sse>& x,
    const simd<double, simd_abi::sse>& y) noexcept {
  return _mm_max_pd(static_cast<__m128d>(x), static_cast<__m128d>(y));
}

} // namespace fub

#endif
#endif // !FUB_SIMD_sse_HPP
