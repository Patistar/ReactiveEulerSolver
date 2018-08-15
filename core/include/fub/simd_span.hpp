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

#include <Vc/vector.h>

namespace fub {
inline namespace v1 {

template <typename T, typename Abi = Vc::VectorAbi::Best<T>,
          std::ptrdiff_t Extent = dynamic_extent,
          typename Alignment = Vc::AlignedTag>
using simd_span = basic_span<T, Extent, accessor_simd<T, Abi, Alignment>>;

template <typename T, typename Abi, typename Extents, typename Layout>
using simd_mdspan =
    basic_mdspan<T, Extents, Layout, accessor_simd<T, Abi, Vc::AlignedTag>>;

} // namespace v1
} // namespace fub

#endif