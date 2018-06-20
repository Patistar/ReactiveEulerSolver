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
#include "fub/patch_view.hpp"
#include "fub/serial/grid.hpp"

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
      dest(indices) = equation.set_momentum(dest(indices), rho_u);
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
      dest(indices) = equation.set_momentum(dest(indices), rho_u);
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
    using traits = grid_traits<Grid>;
    using node_type = typename traits::node_type;
    using patch_type = typename traits::patch_type;
    using variables_tuple = typename patch_type::variables_tuple;
    using equation_type = typename traits::equation_type;
    using extents_type = typename traits::extents_type;
    using result_extents_type =
        replace_extent_t<extents_type, as_int(Axis), int_c<Width>>;
    using result_node_type =
        typename traits::template bind_node<equation_type, result_extents_type>;
    auto where = traits::locality(partition);
    node_type node = traits::node(partition);
    auto reflect_data = [](const node_type& left,
                           const equation_type& equation) {
      result_node_type result = grid_traits<Grid>::dataflow(
          [](auto view, auto where, const equation_type& equation) {
            auto left = view.get();
            result_extents_type reduced_extents = replace_extent(
                left.extents(), int_c<as_int(Axis)>, int_c<Width>);
            auto right = make_patch(variables_tuple(), reduced_extents);
            apply<Axis, Dir>(left, make_view(right), equation);
            return result_node_type(where.get(), std::move(right));
          },
          left.get_patch_view(), left.get_locality(),
          equation);
      return result;
    };
    result_node_type result = grid_traits<Grid>::dataflow_action(
        reflect_data, std::move(where), node, grid.equation());
    return result;
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
} // namespace fub

#endif
