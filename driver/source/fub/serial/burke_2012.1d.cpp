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

#include "fub/serial/burke_2012.1d.hpp"
#include "fub/serial/advance.hpp"
#include "fub/serial/initialise.hpp"

#include "fub/euler/hlle_riemann_solver.hpp"
#include "fub/godunov_method.hpp"
#include "fub/hyperbolic_system_solver.hpp"
#include "fub/patch_view.hpp"
#include "fub/time_integrator/forward_euler.hpp"

namespace fub {
namespace serial {
namespace {
const fub::godunov_method<fub::euler::hlle_riemann_solver> flux_method;

const fub::time_integrator::forward_euler time_integrator;

const fub::hyperbolic_system_solver<
    serial::burke_2012_1d::grid_type, serial::burke_2012_1d::boundary_condition,
    uniform_cartesian_coordinates<1>, decltype(flux_method),
    decltype(time_integrator)>
    advective_solver{flux_method, time_integrator};
} // namespace

burke_2012_1d::state_type burke_2012_1d::initialise(
    initial_condition_function f,
    const uniform_cartesian_coordinates<rank>& coordinates, int depth) {
  return serial::initialise<state_type>(std::move(f), coordinates, depth, 0.25);
}

burke_2012_1d::state_type burke_2012_1d::advance(
    const state_type& state, std::chrono::duration<double> goal,
    const boundary_condition& boundary, feedback_function feedback) {
  return serial::advance(advective_solver, state, goal, boundary,
                         std::move(feedback));
}

} // namespace serial
} // namespace fub
