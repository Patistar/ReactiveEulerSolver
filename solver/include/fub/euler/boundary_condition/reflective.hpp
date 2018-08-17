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

#ifndef FUB_EULER_REFLECTIVE_HPP
#define FUB_EULER_REFLECTIVE_HPP

#include "fub/algorithm.hpp"

#include "fub/equation.hpp"
#include "fub/euler/variables.hpp"
#include "fub/face.hpp"
#include "fub/grid.hpp"

#include <cassert>

namespace fub {
inline namespace v1 {
namespace euler {
namespace boundary_condition {

template <std::size_t N>
std::array<std::ptrdiff_t, N>
reflect(std::array<std::ptrdiff_t, N> index, int dim,
        std::array<std::ptrdiff_t, N> extents) noexcept {
  index[dim] = extents[dim] - index[dim] - 1;
  return index;
}

class reflective {
public:
  template <typename Equation, typename Grid, typename Quadrant>
  void compute_neighbor_data(Equation equation, const Grid& grid,
                             const Quadrant& quad, face f,
                             patch_t<Equation> patch) {
    static constexpr int Rank = Equation::rank();
    auto&& origin = grid.patch_data(quad);
    assert(origin.get_extents() == patch.get_extents());
    const std::array<std::ptrdiff_t, Rank> e = as_array(patch.get_extents());
    for_each_index(origin.get_mapping(), [=](auto... is) {
      const std::array<std::ptrdiff_t, Rank> index{{is...}};
      const std::array<std::ptrdiff_t, Rank> reflected =
          reflect(index, as_int(f.dimension), e);
      patch(index) = origin(reflected);
      auto dyn_mom = dynamic_momentum<Rank>(as_int(f.dimension));
      visit([&](auto mom) { patch[mom](index) *= -1; }, dyn_mom);
    });
  }
};

constexpr bool operator==(const reflective&, const reflective&) noexcept {
  return true;
}
constexpr bool operator!=(const reflective&, const reflective&) noexcept {
  return false;
}

} // namespace boundary_condition
} // namespace euler
} // namespace v1
} // namespace fub

#endif
