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
#include "fub/patch_view.hpp"

#include <cassert>

namespace fub {
namespace euler {
namespace boundary_condition {

class reflective {
  template <axis Axis, direction Dir, typename Source, typename Dest,
            typename Equation>
  static void apply(const Source& source, const Dest& dest,
                    const Equation& equation, std::true_type) noexcept {
    static constexpr int Dim = as_int(Axis);
    for_each_index(dest.extents(), [&](auto indices) {
      const std::ptrdiff_t size = source.extents().get(Dim);
      const std::ptrdiff_t i = indices[Dim];
      auto mapped = indices;
      mapped[Dim] = size - 1 - i;
      dest(indices) = source(mapped);
      auto rho_u = equation.get_momentum(dest(indices));
      rho_u[Dim] *= -1;
      equation.set_momentum(dest(indices), rho_u);
    });
  }

  template <axis Axis, direction Dir, typename Source, typename Dest,
            typename Equation>
  static void apply(const Source& source, const Dest& dest,
                    const Equation& equation, std::false_type) noexcept {
    using namespace variables;
    static constexpr int Dim = as_int(Axis);
    for_each_index(dest.extents(), [&](auto indices) {
      const std::ptrdiff_t size = dest.extents().get(Dim);
      const std::ptrdiff_t i = indices[Dim];
      auto mapped = indices;
      mapped[Dim] = size - 1 - i;
      dest(mapped) = source(indices);
      auto rho_u = equation.get_momentum(dest(indices));
      rho_u[Dim] *= -1;
      equation.set_momentum(dest(indices), rho_u);
    });
  }

  template <axis Axis, direction Dir, typename Source, typename Dest,
            typename Equation>
  static void apply(const Source& source, const Dest& dest,
                    const Equation& equation) noexcept {
    return apply<Axis, Dir>(source, dest, equation,
                            bool_c<Dir == direction::right>);
  }

public:
  template <int Width, axis Axis, direction Dir, typename Grid,
            typename Coordinates>
  static auto
  get_face_neighbor(const typename grid_traits<Grid>::partition_type& partition,
                    const Grid& grid, const Coordinates& /* coordinates */) {
    using B = typename grid_traits<Grid>::patch_type;
    using Variables = typename B::variables_tuple;
    auto reflect_data = [equation = grid.equation()](const B& left) {
      auto reduced_extents =
          replace_extent(left.extents(), int_c<as_int(Axis)>, int_c<Width>);
      auto right = make_patch(Variables(), reduced_extents, left.descriptor());
      apply<Axis, Dir>(make_view(left), make_view(right), equation);
      return right;
    };
    return grid_traits<Grid>::dataflow(reflect_data, partition);
  }
};

const bool operator==(const reflective&, const reflective&) noexcept {
  return true;
}
const bool operator!=(const reflective&, const reflective&) noexcept {
  return false;
}

} // namespace boundary_condition
} // namespace euler
} // namespace fub

#endif
