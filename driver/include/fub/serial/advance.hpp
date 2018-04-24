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

#ifndef FUB_SOLVER_SERIAL_ADVANCE_HPP
#define FUB_SOLVER_SERIAL_ADVANCE_HPP

#include "fub/algorithm.hpp"
#include "fub/functional.hpp"
#include "fub/ode_solver/radau.hpp"
#include <chrono>
#include <utility>

namespace fub {
namespace serial {

template <typename Solver, typename State, typename BoundaryCondition,
          typename Feedback>
State advance(const Solver &solver, State state,
              std::chrono::duration<double> goal,
              BoundaryCondition boundary_condition, Feedback feedback) {
  while (state.time < goal) {
    std::chrono::duration<double> dt = goal - state.time;
    decltype(state.grid) grid;
    std::tie(grid, dt) = solver.get_next_time_step(
        state.grid, state.cfl, dt, state.coordinates, boundary_condition);
    state.grid = std::move(grid);
    state.time = state.time + dt;
    state.dt = dt;
    state.cycle = state.cycle + 1;
    if (!::fub::invoke(feedback, state)) {
      break;
    }
  }
  return state;
}

} // namespace serial
} // namespace fub

#endif // !HPX_ADVANCE_HPP
