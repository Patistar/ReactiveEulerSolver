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

#include "fub/single_stage_1d_serial.hpp"
#include "fub/serial_advance.hpp"
#include "fub/initialise.hpp"

#include "fub/euler/boundary_condition/reflective.hpp"
#include "fub/euler/hlle_riemann_solver.hpp"
#include "fub/euler/muscl_hancock_method.hpp"
#include "fub/hyperbolic_system_solver.hpp"
#include "fub/patch_view.hpp"
#include "fub/time_integrator/forward_euler.hpp"

namespace fub {
namespace {
const single_stage_1d_serial::equation_type equation{};

const fub::euler::muscl_hancock_method<fub::euler::hlle_riemann_solver>
    flux_method;

const fub::time_integrator::forward_euler time_integrator;

const fub::hyperbolic_system_solver<decltype(equation), decltype(flux_method),
                                    decltype(time_integrator)>
    advective_solver{equation, flux_method, time_integrator};

const fub::euler::boundary_condition::reflective boundary_condition;
} // namespace

single_stage_1d_serial::state_type single_stage_1d_serial::initialise(
    initial_condition_function f,
    const uniform_cartesian_coordinates<rank>& coordinates, int depth) {
  using traits = grid_traits<grid_type>;
  grid_type grid(depth, patch_extents_type());
  for (auto& partition : grid) {
    auto octant = traits::octant(partition);
    partition.second = traits::dataflow(
        [=, g = std::ref(f)](patch_type patch) {
          auto adapted = adapt(coordinates, octant);
          auto view = make_view(patch);
          for_each_index(patch.extents(),
                         [&](const std::array<index, rank>& i) {
                           view(i) = fub::invoke(g, fub::apply(adapted, i));
                         });
          return std::make_shared<traits::node_type>(std::move(patch));
        },
        std::move(partition));
  }
  return state_type{std::move(grid),
                    coordinates,
                    std::chrono::duration<double>(0),
                    std::chrono::duration<double>(0),
                    0,
                    0.8};
}

single_stage_1d_serial::state_type
single_stage_1d_serial::advance(const state_type& state,
                                std::chrono::duration<double> goal,
                                feedback_function feedback) {
  return serial_advance(advective_solver, state, goal, boundary_condition,
                        std::move(feedback));
}

} // namespace fub
