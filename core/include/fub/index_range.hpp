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

/// @file This file introduces `extents<E0, ..., En>`, a compact
/// multi-dimensional size type. Each integral extent `Ei` stands an upper bound
/// in dimension `i` and can either be a compile-time constant signed integral
/// value or `dynamic_extent`. Each compile-time sized extent does not take any
/// extra byte.

#ifndef FUB_CORE_INDEX_RANGE_HPP
#define FUB_CORE_INDEX_RANGE_HPP

#include "fub/extents.hpp"
#include "fub/layout_left.hpp"
#include "fub/type_traits.hpp"

#include <array>
#include <cassert>

namespace fub {
inline namespace v1 {

template <int Rank> struct index_range {
  std::array<index, Rank> origin;
  std::array<index, Rank> extents;
};

template <int Rank, typename F>
void for_each_index(index_range<Rank> range, F feedback) {
  dynamic_extents_t<Rank> extents(range.extents);
  using mapping_t =
      typename layout_left::template mapping<dynamic_extents_t<Rank>>;
  for_each_index(mapping_t(extents), [=](std::array<index, Rank> indices) {
    std::transform(indices.begin(), indices.end(), range.origin.begin(),
                   indices.begin(), [=](index i, index o) { return i - o; });
    fub::invoke(feedback, indices);
  });
}

} // namespace v1
} // namespace fub

#endif