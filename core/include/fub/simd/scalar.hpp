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

#ifndef FUB_SIMD_SCALAR_HPP
#define FUB_SIMD_SCALAR_HPP

#include "fub/simd/simd_fwd.hpp"

namespace fub {

////////////////////////////////////////////////////////////////////////////////
//                                                    Scalar SIMD Specialisation

/// @brief Simd wrapper for the scalar case (SISD case).
template <typename T> class simd<T, simd_abi::scalar> {
public:
  using value_type = T;
  using abi_type = simd_abi::scalar;
  using mask_type = simd_mask<value_type, abi_type>;
  using native_type = T;

  // Constructors

  simd() = default;

  /// Implicit conversion from value.
  /// 
  /// Effect: Initialises m_value with u.
  ///
  /// Throws: nothing
  ///
  /// Postcondition: m_value == u
  template <
      typename U,
      std::enable_if_t<std::is_convertible<U, value_type>::value>* = nullptr>
  simd(U&& u) noexcept : m_value(u) {}

  /// Construct by load from memory location.
  /// 
  /// Effect: As if
  ///     simd<T, simd_abi::scalar> v;
  ///     v.copy_from(memory, flag)
  ///
  /// Throws: nothing
  ///
  /// Postcondition: m_value == *memory
  template <
      typename U, typename AlignmentFlag,
      std::enable_if_t<std::is_convertible<U, value_type>::value>* = nullptr,
      std::enable_if_t<is_simd_flag_type<AlignmentFlag>::value>* = nullptr>
  simd(const U* memory, AlignmentFlag /* flag */) noexcept : m_value(*memory) {}

  // Observers

  static constexpr int size() noexcept { return 1; }

  /// Load a value from memory.
  ///
  /// Note: There are no alignment restrictions in the scalar case.
  ///
  /// Throws: Nothing.
  template <
      typename U, typename AlignmentFlag,
      std::enable_if_t<std::is_convertible<U, value_type>::value>* = nullptr,
      std::enable_if_t<is_simd_flag_type<AlignmentFlag>::value>* = nullptr>
  constexpr void
  copy_from(const U* mem, Flags) noexcept {
    m_value = *mem;
  }

  /// Store a value to memory.
  ///
  /// Note: There are no alignment restrictions in the scalar case.
  ///
  /// Throws: Nothing.
  template <typename U, typename Flags>
  constexpr std::enable_if_t<std::is_convertible<value_type, U>::value>
  copy_to(U* mem, Flags) const noexcept {
    *mem = static_cast<U>(m_value);
  }

  // Element Access
#if __cplusplus >= 201703L
  constexpr value_type operator[]([[maybe_unused]] int i) const noexcept {
#elif defined(NDEBUG)
  constexpr value_type operator[](int) const noexcept {
#else
  constexpr value_type operator[](int i) const noexcept {
#endif
    assert(i == 0);
    return m_value;
  }

  // @brief This makes the simd-type static_cast-able to its underlying native
  // handle.
  explicit constexpr operator native_type() const noexcept { return m_value; }

  // Arithmetic Operators

  friend constexpr simd operator+(const simd& x, const simd& y) noexcept {
    return x.m_value + y.m_value;
  }
  friend constexpr simd operator-(const simd& x, const simd& y) noexcept {
    return x.m_value - y.m_value;
  }
  friend constexpr simd operator/(const simd& x, const simd& y) noexcept {
    return x.m_value / y.m_value;
  }
  friend constexpr simd operator*(const simd& x, const simd& y) noexcept {
    return x.m_value * y.m_value;
  }

  constexpr simd operator-() const noexcept {
    return simd{-m_value};
  }

  constexpr simd& operator+=(const simd& x) noexcept {
    m_value += x.m_value;
    return *this;
  }
  constexpr simd& operator-=(const simd& x) noexcept {
    m_value -= x.m_value;
    return *this;
  }

  // Comparison Operators

  friend mask_type operator<(const simd& x, const simd& y) noexcept {
    return mask_type{x.m_value < y.m_value};
  }
  friend mask_type operator<=(const simd& x, const simd& y) noexcept {
    return mask_type{x.m_value <= y.m_value};
  }
  friend mask_type operator>(const simd& x, const simd& y) noexcept {
    return mask_type{x.m_value > y.m_value};
  }
  friend mask_type operator>=(const simd& x, const simd& y) noexcept {
    return mask_type{x.m_value >= y.m_value};
  }

  friend mask_type operator==(const simd& x, const simd& y) noexcept {
    return mask_type{x.m_value == y.m_value};
  }
  friend mask_type operator!=(const simd& x, const simd& y) noexcept {
    return mask_type{x.m_value != y.m_value};
  }

private:
  native_type m_value{};
};

template <typename T> class simd_mask<T, simd_abi::scalar> {
public:
  using value_type = bool;
  using abi_type = simd_abi::scalar;
  using simd_type = simd<T, abi_type>;

  simd_mask() = default;
  constexpr explicit simd_mask(value_type value) : m_value{value} {}
  explicit constexpr operator value_type() const noexcept { return m_value; }

  // Load and Store from Memory

  template <typename Flags>
  constexpr std::enable_if_t<is_simd_flag_type_v<Flags>, void>
  copy_from(const value_type* mem, Flags) noexcept {
    m_value = *mem;
  }

  template <typename Flags>
  constexpr std::enable_if_t<is_simd_flag_type_v<Flags>, void>
  copy_to(value_type* mem, Flags) noexcept {
    *mem = m_value;
  }

  // Element Access

#if __cplusplus >= 201703L
  constexpr value_type operator[]([[maybe_unused]] int i) const noexcept {
#elif defined(NDEBUG)
  constexpr value_type operator[](int) const noexcept {
#else
  constexpr value_type operator[](int i) const noexcept {
#endif
    assert(i == 0);
    return m_value;
  }

  // Unary Operators (element-wise)

  constexpr simd_mask operator!() const noexcept { return simd_mask{!m_value}; }

  // Binary Operators (element-wise)

  friend constexpr simd_mask operator&&(const simd_mask& m1,
                                        const simd_mask& m2) noexcept {
    return simd_mask{m1.m_value && m2.m_value};
  }

  friend constexpr simd_mask operator||(const simd_mask& m1,
                                        const simd_mask& m2) noexcept {
    return simd_mask{m1.m_value || m2.m_value};
  }

  friend constexpr simd_mask operator&(const simd_mask& m1,
                                       const simd_mask& m2) noexcept {
    return simd_mask{m1.m_value & m2.m_value};
  }

  friend constexpr simd_mask operator|(const simd_mask& m1,
                                       const simd_mask& m2) noexcept {
    return simd_mask{m1.m_value | m2.m_value};
  }

  friend constexpr simd_mask operator^(const simd_mask& m1,
                                       const simd_mask& m2) noexcept {
    return simd_mask{m1.m_value ^ m2.m_value};
  }

  // Compound Operators (element-wise)

  friend constexpr simd_mask& operator&=(simd_mask& m1,
                                         const simd_mask& m2) noexcept {
    m1.value &= m2.m_value;
    return m1;
  }

  friend constexpr simd_mask& operator|=(simd_mask& m1,
                                         const simd_mask& m2) noexcept {
    m1.m_value |= m2.m_value;
    return m1;
  }

  friend constexpr simd_mask& operator^=(simd_mask& m1,
                                         const simd_mask& m2) noexcept {
    m1.m_value ^= m2.m_value;
    return m1;
  }

  // Compare Operators (element-wise)

  friend constexpr simd_mask operator==(const simd_mask& m1,
                                        const simd_mask& m2) noexcept {
    return simd_mask{m1.m_value == m2.m_value};
  }

  friend constexpr simd_mask operator!=(const simd_mask& m1,
                                        const simd_mask& m2) noexcept {
    return simd_mask{m1.m_value != m2.m_value};
  }

private:
  bool m_value;
};

template <typename T> using scalar = simd<T, simd_abi::scalar>;

template <typename T>
bool any_of(const simd_mask<T, simd_abi::scalar>& x) noexcept {
  return static_cast<bool>(x);
}

template <typename T>
bool all_of(const simd_mask<T, simd_abi::scalar>& x) noexcept {
  return static_cast<bool>(x);
}

template <typename T>
class const_where_expression<simd_mask<T, simd_abi::scalar>,
                             simd<T, simd_abi::scalar>> {
public:
  using value_type = T;
  using abi_type = simd_abi::scalar;
  using mask_type = simd_mask<T, simd_abi::scalar>;
  using simd_type = simd<T, simd_abi::scalar>;

  const_where_expression(const const_where_expression&);
  const_where_expression& operator=(const const_where_expression&);

  simd_type operator-() const&& { return m_mask ? -m_simd : m_simd; }

  template <typename U, typename Flags>
  std::enable_if_t<
      is_simd_flag_type_v<Flags> && std::is_convertible<T, U>::value, void>
  copy_to(U* mem, Flags) const&& noexcept {
    if (m_mask) {
      m_simd.copy_to(mem);
    }
  }

private:
  const mask_type& m_mask{};
  simd_type& m_simd{};

  const_where_expression(const mask_type& m, simd_type& s) noexcept
      : m_mask{m}, m_simd{s} {}

  friend class where_expression<mask_type, simd_type>;

public:
  friend const_where_expression where(const mask_type& mask,
                                      const simd_type& simd) noexcept {
    return const_where_expression(mask, const_cast<simd_type&>(simd));
  }
};

template <typename T>
class where_expression<simd_mask<T, simd_abi::scalar>,
                       simd<T, simd_abi::scalar>>
    : public const_where_expression<simd_mask<T, simd_abi::scalar>,
                                    simd<T, simd_abi::scalar>> {
public:
  using value_type = T;
  using abi_type = simd_abi::scalar;
  using mask_type = simd_mask<T, simd_abi::scalar>;
  using simd_type = simd<T, simd_abi::scalar>;

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
      is_simd_flag_type_v<Flags> && std::is_convertible<U, T>::value, void>
  copy_from(const U* mem, Flags) const&& noexcept {
    if (m_mask) {
      m_simd.copy_from(mem);
    }
  }

  template <typename U>
  std::enable_if_t<std::is_constructible<simd_type, U>::value, void>
  operator=(U&& u) {
    if (m_mask) {
      m_simd = std::forward<U>(u);
    }
  }
};

////////////////////////////////////////////////////////////////////////////////
// Math Function Overloads

template <typename T>
simd<T, simd_abi::scalar> sqrt(const simd<T, simd_abi::scalar>& x) noexcept {
  return std::sqrt(static_cast<T>(x));
}

template <typename T>
simd<T, simd_abi::scalar> abs(const simd<T, simd_abi::scalar>& x) noexcept {
  return x < 0 ? simd<T, simd_abi::scalar>(-static_cast<T>(x)) : x;
}

template <typename T>
simd<T, simd_abi::scalar> exp(const simd<T, simd_abi::scalar>& x) noexcept {
  return std::exp(static_cast<T>(x));
}

template <typename T>
simd<T, simd_abi::scalar> log(const simd<T, simd_abi::scalar>& x) noexcept {
  return std::log(static_cast<T>(x));
}

template <typename T>
simd<T, simd_abi::scalar> round(const simd<T, simd_abi::scalar>& x) noexcept {
  return std::round(static_cast<T>(x));
}

////////////////////////////////////////////////////////////////////////////////
// Special Math Function Overloads

template <typename T>
simd<T, simd_abi::scalar> min(const simd<T, simd_abi::scalar>& x,
                              const simd<T, simd_abi::scalar>& y) noexcept {
  return std::min(static_cast<T>(x), static_cast<T>(y));
}

template <typename T>
simd<T, simd_abi::scalar> max(const simd<T, simd_abi::scalar>& x,
                              const simd<T, simd_abi::scalar>& y) noexcept {
  return std::max(static_cast<T>(x), static_cast<T>(y));
}
} // namespace fub

#endif // !FUB_SIMD_SCALAR_HPP
