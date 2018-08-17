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

#ifndef FUB_SOLVER_PERMUTATE_DIMENSIONS_HPP
#define FUB_SOLVER_PERMUTATE_DIMENSIONS_HPP

#include "fub/variable_view.hpp"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#include <unsupported/Eigen/CXX11/Tensor>
#pragma clang diagnostic pop

#include <algorithm>

namespace fub {
inline namespace v1 {
template <std::size_t Rank>
std::array<int, Rank> _permutation_index(std::array<int, 2> p) noexcept {
  std::array<int, Rank> is;
  std::iota(is.begin(), is.end(), 0);
  std::swap(is[p[0]], is[p[1]]);
  return is;
}

template <typename VL, typename MdSpan>
void permutate_dimensions(
    basic_variable_view<VL, MdSpan> in, const std::array<int, 2>& permutation,
    nodeduce_t<basic_variable_view<VL, mdspan_remove_const_t<MdSpan>>>
        out) noexcept {
  constexpr int Rank = MdSpan::rank();
  auto extents = as_array(in.get_extents());
  using T = typename MdSpan::value_type;
  using Tensor = Eigen::TensorMap<Eigen::Tensor<T, Rank>>;
  using ConstTensor = Eigen::TensorMap<Eigen::Tensor<const T, Rank>>;
  for_each(in.get_variable_list(), [&](auto var) {
    ConstTensor view_in(in[var].span().data(), extents);
    Tensor view_out(out[var].span().data(), extents);
    view_out = view_in.shuffle(_permutation_index<Rank>(permutation));
  });
}

} // namespace v1
} // namespace fub

#endif