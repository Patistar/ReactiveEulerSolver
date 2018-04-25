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

#include "fub/hpx/kinetic.gri_30.1d.hpp"
#include "fub/hpx/advance.hpp"
#include "fub/hpx/initialise.hpp"
#include "fub/patch_view.hpp"

#include "fub/euler/boundary_condition/reflective.hpp"
#include "fub/euler/hlle_riemann_solver.hpp"
#include "fub/euler/kinetic_source_term.hpp"
#include "fub/euler/muscl_hancock_method.hpp"
#include "fub/strang_splitting.hpp"
#include "fub/hyperbolic_system_solver.hpp"
#include "fub/hyperbolic_system_source_solver.hpp"
#include "fub/patch_view.hpp"
#include "fub/time_integrator/forward_euler.hpp"

namespace fub {
namespace hpx {
namespace kinetic {
namespace {
const gri_30_1d::equation_type equation{};
const euler::muscl_hancock_method<euler::hlle_riemann_solver> flux_method;
const time_integrator::forward_euler time_integrator;
const hyperbolic_system_solver<decltype(equation), decltype(flux_method),
                               decltype(time_integrator)>
    advective_solver{equation, flux_method, time_integrator};
const auto kinetic_source_term = euler::make_kinetic_source_term(equation);
const auto solver = make_hyperbolic_system_source_solver(
    strang_splitting(), advective_solver, kinetic_source_term);
} // namespace

gri_30_1d::state_type
gri_30_1d::initialise(initial_condition_function f,
                      const uniform_cartesian_coordinates<rank>& coordinates,
                      int depth) {
  return ::fub::hpx::initialise<state_type>(std::move(f), coordinates, depth);
}

gri_30_1d::state_type gri_30_1d::advance(const state_type& state,
                                         std::chrono::duration<double> goal,
                            const boundary_condition& boundary,
                                         feedback_function feedback) {
  return ::fub::hpx::advance(solver, state, goal, boundary,
                             std::move(feedback));
}

} // namespace kinetic
} // namespace hpx
} // namespace fub
