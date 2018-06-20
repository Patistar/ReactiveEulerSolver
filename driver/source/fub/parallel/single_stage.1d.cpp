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

#include "fub/parallel/single_stage.1d.hpp"
#include "fub/parallel/advance.hpp"
#include "fub/parallel/initialise.hpp"
#include "fub/patch_view.hpp"

#include "fub/euler/boundary_condition/reflective.hpp"
#include "fub/euler/hlle_riemann_solver.hpp"
#include "fub/euler/muscl_hancock_method.hpp"
#include "fub/hyperbolic_system_solver.hpp"
#include "fub/time_integrator/forward_euler.hpp"

namespace fub {
namespace parallel {
namespace {
const godunov_method<fub::euler::hlle_riemann_solver> flux_method;

const fub::time_integrator::forward_euler time_integrator;

const fub::hyperbolic_system_solver<parallel::single_stage_1d::grid_type,
                                    parallel::single_stage_1d::boundary_condition,
                                    uniform_cartesian_coordinates<1>,
                                    decltype(flux_method),
                                    decltype(time_integrator)>
    advective_solver{flux_method, time_integrator};
} // namespace

single_stage_1d::state_type single_stage_1d::initialise(
    initial_condition_function f,
    const uniform_cartesian_coordinates<rank>& coordinates, int depth) {
  return parallel::initialise<state_type>(std::move(f), coordinates, depth);
}

single_stage_1d::state_type single_stage_1d::advance(
    const state_type& state, std::chrono::duration<double> goal,
    const boundary_condition& boundary, feedback_function feedback) {
  return parallel::advance(advective_solver, state, goal, boundary,
                             std::move(feedback));
}

} // namespace parallel
} // namespace fub
