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

#ifndef FUB_CORE_SIMD_SPAN_HPP
#define FUB_CORE_SIMD_SPAN_HPP

#include "fub/accessor_simd.hpp"
#include "fub/layout_simd.hpp"
#include "fub/mdspan.hpp"
#include "fub/simd.hpp"

namespace fub {
namespace simd_span_detail {
/// Returns by default a simd span with an unaligned accessor.
template <typename T, typename Abi, std::ptrdiff_t Extent, typename Alignment>
struct simd_span_detector {
  using type = basic_mdspan<T, extents<Extent>, layout_simd_padded<T, Abi>,
                            accessor_simd_unaligned<T, Abi>>;
};

/// This specialization returns a simd span with an aligned accessor if
/// `vector_alignment_tag` is specified.
template <typename T, typename Abi, std::ptrdiff_t Extent>
struct simd_span_detector<T, Abi, Extent, flags::vector_aligned_tag> {
  using type = basic_mdspan<T, extents<Extent>, layout_simd_padded<T, Abi>,
                            accessor_simd_aligned<T, Abi>>;
};
} // namespace simd_span_detail

template <typename T, typename Abi = simd_abi::native<remove_cvref_t<T>>,
          std::ptrdiff_t Extent = dynamic_extent,
          typename Alignment = flags::vector_aligned_tag>
using simd_span =
    typename simd_span_detail::simd_span_detector<T, Abi, Extent,
                                                  Alignment>::type;

} // namespace fub

#endif