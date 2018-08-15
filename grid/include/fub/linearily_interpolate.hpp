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

#ifndef FUB_GRID_LINEARILY_INTERPOLATE_HPP
#define FUB_GRID_LINEARILY_INTERPOLATE_HPP

#include "fub/index_range.hpp"
#include "fub/p4est/quadrant.hpp"
#include "fub/variable_view.hpp"

namespace fub {
inline namespace v1 {
constexpr std::array<int, 2> relative_coordinates(int child, int_constant<2>) {
  return {child % 2, child / 2};
}

constexpr std::array<int, 3> relative_coordinates(int child, int_constant<3>) {
  return {child % 2, child / 2, child / 4};
}

template <std::size_t Rank>
index_range<Rank> sub_index_range(std::array<index, Rank> extents, int child) {
  std::transform(extents.begin(), extents.end(), extents.begin(),
                 [](index i) { return i / 2; });
  std::array<index, Rank> coords = relative_coordinates(child, int_c<Rank>);
  std::array<index, Rank> origin;
  std::transform(extents.begin(), extents.end(), coords.begin(), origin.begin(),
                 [](index extent, index coord) { return extent * coord; });
  return {origin, extents};
}

template <typename VL, typename T, typename E, typename A>
void linearily_refine(variable_view<VL, const T, E, A> coarse,
                      nodeduce_t<variable_view<VL, T, E, rebind_t<A, T>>> fine,
                      int child) noexcept {
  constexpr int Rank = Coarse::rank();
  const index_range<Rank> box = sub_index_range(coarse.get_extents(), child);
  for_each_index(box, [=](std::array<index, Rank> index) {
    const auto refined = refine(index - box.origin);
    const auto data = coarse(index);
    std::for_each(refined.begin(), refined.end(),
                  [=](std::array<index, Rank> r) { fine(r) = data; });
  });
}

template <typename VL, typename T, typename E, typename A>
void linearily_coarsen(
    variable_view<VL, const T, E, A> fine,
    nodeduce_t<variable_view<VL, T, E, rebind_t<A, T>>> coarse,
    int child) noexcept {
  constexpr int Rank = Fine::rank();
  const index_range<Rank> box = sub_index_range(coarse.get_extents(), child);
  for_each_index(box, [=](std::array<index, Rank> index) {
    const auto refined = refine(index - box.origin);
    coarse(index) = fine(refined[0]);
    std::for_each(std::next(refined.begin()), refined.end(),
                  [=](std::array<index, Rank> r) {
                    fub::for_each(coarse.get_variable_list(), [=](auto var) {
                      coarse[var](index) += fine[var](r);
                    });
                  });
    fub::for_each(coarse.get_variable_list(),
                  [=](auto var) { coarse[var](index) /= 4; });
  });
}

} // namespace v1
} // namespace fub

#endif