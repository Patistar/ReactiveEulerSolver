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

#ifndef FUB_SOLVER_MAIN_DRIVER_HPP
#define FUB_SOLVER_MAIN_DRIVER_HPP

#include "fub/algorithm.hpp"
#include "fub/functional.hpp"
#include "fub/print_cycle_timings.hpp"

#include <cmath>
#include <stdexcept>

namespace fub {

struct simulation_options {
  std::chrono::duration<double> feedback_interval;
  std::chrono::duration<double> final_time;
};

template <typename Solver, typename BoundaryCondition,
          typename IntervalFeedback,
          typename CycleFeedback = print_cycle_timings>
void run_simulation(
    const Solver& solver, typename Solver::state_type state,
    const BoundaryCondition& boundary_condition,
    const simulation_options& options,
    const IntervalFeedback& interval_feedback,
    const CycleFeedback& cycle_feedback = CycleFeedback()) noexcept {
  constexpr std::chrono::duration<double> eps(
      std::numeric_limits<double>::epsilon());
  while (state.time < options.final_time) {
    const std::chrono::duration<double> rest =
        options.final_time - state.time + eps;
    const std::chrono::duration<double> dt =
        std::min(options.feedback_interval, rest);
    state = solver.advance(std::move(state), state.time + dt,
                           boundary_condition, std::ref(cycle_feedback));
    fub::invoke(interval_feedback, state);
  }
}

} // namespace fub

#endif
