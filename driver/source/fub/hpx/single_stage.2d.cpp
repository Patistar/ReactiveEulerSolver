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

#include "fub/hpx/single_stage.2d.hpp"
#include "fub/hpx/advance.hpp"
#include "fub/hpx/initialise.hpp"
#include "fub/patch_view.hpp"

#include "fub/euler/boundary_condition/reflective.hpp"
#include "fub/euler/hlle_riemann_solver.hpp"
#include "fub/euler/muscl_hancock_method.hpp"
#include "fub/hyperbolic_system_solver.hpp"
#include "fub/time_integrator/forward_euler.hpp"

namespace fub {
namespace hpx {
namespace {
const single_stage_2d::equation_type equation{};

const fub::euler::muscl_hancock_method<fub::euler::hlle_riemann_solver>
    flux_method;

const fub::time_integrator::forward_euler time_integrator;

const fub::hyperbolic_system_solver<decltype(equation), decltype(flux_method),
                                    decltype(time_integrator)>
    advective_solver{equation, flux_method, time_integrator};
} // namespace

single_stage_2d::state_type single_stage_2d::initialise(
    initial_condition_function f,
    const uniform_cartesian_coordinates<rank>& coordinates, int depth) {
  return ::fub::hpx::initialise<state_type>(std::move(f), coordinates, depth);
}

namespace {
template <class BX, class BY> struct joined_boundary_condition {
  BX m_boundary_x;
  BY m_boundary_y;

  template <int Width, direction Dir, typename Grid, typename Coordinates>
  auto get_face_neighbor_impl(
      const typename grid_traits<Grid>::partition_type& partition,
      const Grid& grid, const Coordinates& coordinates, int_constant<0>) const {
    return m_boundary_x.template get_face_neighbor<Width, axis::x, Dir>(
        partition, grid, coordinates);
  }

  template <int Width, direction Dir, typename Grid, typename Coordinates>
  auto get_face_neighbor_impl(
      const typename grid_traits<Grid>::partition_type& partition,
      const Grid& grid, const Coordinates& coordinates, int_constant<1>) const {
    return m_boundary_y.template get_face_neighbor<Width, axis::y, Dir>(
        partition, grid, coordinates);
  }

  template <int Width, axis Axis, direction Dir, typename Grid,
            typename Coordinates>
  auto
  get_face_neighbor(const typename grid_traits<Grid>::partition_type& partition,
                    const Grid& grid, const Coordinates& coordinates) const {
    return get_face_neighbor_impl<Width, Dir>(partition, grid, coordinates,
                                              int_c<as_int(Axis)>);
  }
};

template <class BX, class BY> auto join_boundary(BX&& bx, BY&& by) {
  return joined_boundary_condition<std::decay_t<BX>, std::decay_t<BY>>{
      std::forward<BX>(bx), std::forward<BY>(by)};
}

} // namespace

single_stage_2d::state_type single_stage_2d::advance(
    const state_type& state, std::chrono::duration<double> goal,
    const boundary_condition_x& boundary_x,
    const boundary_condition_y& boundary_y, feedback_function feedback) {
  return ::fub::hpx::advance(advective_solver, state, goal,
                             join_boundary(boundary_x, boundary_y),
                             std::move(feedback));
}

} // namespace hpx
} // namespace fub
