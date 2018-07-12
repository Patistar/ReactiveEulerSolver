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

#ifndef FUB_SIMD_SIMD_FWD_HPP
#define FUB_SIMD_SIMD_FWD_HPP

#include <array>
#include <cassert>
#include <cmath>
#include <type_traits>

namespace fub {
////////////////////////////////////////////////////////////////////////////////
//                                  Forward Declerations for SIMD related types.

namespace simd_abi {
struct scalar {};
template <typename T> struct native_abi { using type = scalar; };
template <typename T> using native = typename native_abi<T>::type;
template <typename T> using compatible = typename native_abi<T>::type;
}; // namespace simd_abi

template <typename T, typename Abi> class simd_mask {
public:
  simd_mask() = delete;
  simd_mask(const simd_mask&) = delete;
  ~simd_mask() = delete;
  simd_mask& operator=(const simd_mask&) = delete;
};

template <typename T, typename Abi = simd_abi::native<T>> class simd {
public:
  using value_type = T;
  using abi_type = Abi;
  using mask_type = simd_mask<T, Abi>;

  simd() = delete;
  simd(const simd&) = delete;
  ~simd() = delete;
  simd& operator=(const simd&) = delete;
};

template <typename M, typename T> class const_where_expression;
template <typename M, typename T> class where_expression;

template <typename T, typename Abi>
where_expression<typename simd<T, Abi>::mask_type, simd<T, Abi>>
where(const typename simd<T, Abi>::mask_type& mask,
      simd<T, Abi>& simd) noexcept {
  return {mask, simd};
}

template <typename T, typename Abi>
const_where_expression<typename simd<T, Abi>::mask_type, const simd<T, Abi>>
where(const typename simd<T, Abi>::mask_type& mask,
      const simd<T, Abi>& simd) noexcept {
  return {mask, simd};
}

struct element_alignment_tag {};
static constexpr element_alignment_tag element_alignment{};

struct vector_alignment_tag {};
static constexpr vector_alignment_tag vector_alignment{};

////////////////////////////////////////////////////////////////////////////////
//                                                            [simd.type_traits]

////////////////////////////////////////////////////////////////////////////////
//                                                    [simd.type_traits.is_simd]

/// Returns true if `std::is_same<Simd, simd<T, Abi>>` for some T and Abi.
/// @{
template <typename Simd> struct is_simd : std::false_type {};

template <typename T, typename Abi>
struct is_simd<simd<T, Abi>> : std::true_type {};

template <typename Simd> static constexpr bool is_simd_v = is_simd<Simd>::value;
/// @}

////////////////////////////////////////////////////////////////////////////////
//                                               [simd.type_traits.is_simd_mask]

template <typename SimdMask> struct is_simd_mask : std::false_type {};
template <typename T, typename Abi>
struct is_simd_mask<simd_mask<T, Abi>> : std::true_type {};
template <typename Simd>
static constexpr bool is_simd_mask_v = is_simd_mask<Simd>::value;

////////////////////////////////////////////////////////////////////////////////
//                                                [simd.type_traits.is_simd_abi]

template <typename Abi> struct is_simd_abi : std::false_type {};
template <> struct is_simd_abi<simd_abi::scalar> : std::true_type {};
template <typename Simd>
static constexpr bool is_simd_abi_v = is_simd_abi<Simd>::value;

////////////////////////////////////////////////////////////////////////////////
//                                          [simd.type_traits.is_simd_flag_type]

template <typename T> struct is_simd_flag_type : std::false_type {};
template <> struct is_simd_flag_type<element_alignment_tag> : std::true_type {};
template <> struct is_simd_flag_type<vector_alignment_tag> : std::true_type {};
template <typename Simd>
static constexpr bool is_simd_flag_type_v = is_simd_flag_type<Simd>::value;

////////////////////////////////////////////////////////////////////////////////
//                                                  [simd.type_traits.simd_size]

template <typename Simd> struct simd_size;
template <typename T, typename Abi>
struct simd_size<simd<T, Abi>>
    : std::integral_constant<int, simd<T, Abi>::size()> {};
template <typename Simd>
static constexpr int simd_size_v = simd_size<Simd>::value;

////////////////////////////////////////////////////////////////////////////////
//                                           [simd.type_traits.memory_alignment]

template <typename Simd, typename U> struct memory_alignment {};
template <typename Simd, typename U>
static constexpr int memory_alignment_v = memory_alignment<Simd, U>::value;

template <typename T, typename U>
struct memory_alignment<simd<T, simd_abi::scalar>, U>
    : std::integral_constant<int, alignof(T)> {};

} // namespace fub

#endif