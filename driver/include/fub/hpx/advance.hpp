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

#ifndef FUB_SOLVER_HPX_ADVANCE_HPP
#define FUB_SOLVER_HPX_ADVANCE_HPP

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#include <hpx/hpx.hpp>
#pragma clang diagnostic pop

#include "fub/functional.hpp"

#include <chrono>
#include <utility>

namespace fub {
namespace hpx {

template <typename Solver, typename State, typename BoundaryCondition,
          typename Feedback>
State advance(const Solver& solver, const State& state,
              std::chrono::duration<double> goal,
              BoundaryCondition boundary_condition, Feedback feedback) {
  assert(state.time < goal);
  using grid_type = decltype(state.grid);
  static constexpr int computation_graph_depth = 5;
  std::array<::hpx::shared_future<State>, computation_graph_depth> time_levels;
  std::array<::hpx::future<std::chrono::duration<double>>,
             computation_graph_depth>
      times;
  auto get_next_level = [=](::hpx::shared_future<State> future_state) {
    State state = future_state.get();
    if (state.time >= goal) {
      return ::hpx::make_ready_future(std::move(state));
    }
    std::chrono::duration<double> max_step_size = goal - state.time;
    ::hpx::future<std::pair<grid_type, std::chrono::duration<double>>> result =
        solver.get_next_time_step(state.grid, state.cfl, max_step_size,
                                  state.coordinates, boundary_condition);
    return result.then([coordinates = state.coordinates, cfl = state.cfl,
                        cycle = state.cycle,
                        old_time = state.time](auto future_result) {
      auto r = future_result.get();
      std::chrono::duration<double> time = old_time + r.second;
      return State{std::move(r.first), coordinates, time,
                   r.second,           cycle + 1,   cfl};
    });
  };
  auto get_time = [](const ::hpx::shared_future<State>& state) {
    return state.get().time;
  };
  time_levels[0] = ::hpx::make_ready_future(state).then(get_next_level);
  times[0] = time_levels[0].then(get_time);
  for (int i = 1; i < computation_graph_depth; ++i) {
    time_levels[i] = time_levels[i - 1].then(get_next_level);
    times[i] = time_levels[i].then(get_time);
  }
  int prev_i = computation_graph_depth - 1;
  int next_i = 0;
  while (true) {
    times[next_i].wait();
    if (!fub::invoke(std::ref(feedback), time_levels[next_i].get())) {
      break;
    }
    time_levels[next_i] = time_levels[prev_i].then(get_next_level);
    times[next_i] = time_levels[next_i].then(get_time);
    prev_i = (prev_i + 1) % computation_graph_depth;
    next_i = (next_i + 1) % computation_graph_depth;
  }
  return time_levels[next_i].get();
}
} // namespace hpx
} // namespace fub

#endif // !HPX_ADVANCE_HPP
