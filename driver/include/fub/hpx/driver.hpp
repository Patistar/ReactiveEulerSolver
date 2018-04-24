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

#ifndef FUB_SOLVER_DRIVER_HPP
#define FUB_SOLVER_DRIVER_HPP

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#include <hpx/hpx.hpp>
#pragma clang diagnostic pop

#include "fub/algorithm.hpp"

#include <boost/program_options.hpp>
#include <fmt/format.h>

#include <chrono>
#include <cmath>
#include <stdexcept>

namespace fub {
namespace hpx {

template <typename Solver, typename Feedback>
void main_driver(boost::program_options::variables_map &vm,
                 const Solver &solver, typename Solver::state_type state,
                 const Feedback &feedback) noexcept {
  try {
    using state_type = typename Solver::state_type;
    // fetch program options
    const int cycles = vm["cycles"].as<int>();
    const int depth = vm["depth"].as<int>();
    const double time = vm["time"].as<double>();
    const double feedback_interval = vm["feedback_interval"].as<double>();
    ::hpx::future<void> write_job = ::hpx::make_ready_future();
    fmt::print("Program Options:\n\tdepth: {}\n\tcycles: {}\n\ttime: "
               "{}s\n\tfeedback_interval: {}s\n",
               depth, cycles, time, feedback_interval);
    const auto wall_start = std::chrono::steady_clock::now();
    while (state.time.count() < time && state.cycle < cycles) {
      std::chrono::duration<double> dt(
          std::min(feedback_interval, time - state.time.count()));
      std::chrono::duration<double> sub_goal = state.time + dt;
      auto start = std::chrono::steady_clock::now();
      auto feedback = [&](const state_type &s) -> bool {
        auto end = std::chrono::steady_clock::now();
        auto wall_time =
            std::chrono::duration_cast<std::chrono::duration<double>>(
                end - wall_start);
        auto diff = std::chrono::duration_cast<std::chrono::duration<double>>(
            end - start);
        fmt::print("[{:.2f}%] T = {:.10f}s, dt = {:.4e}s, cycle = {:10}   "
                   "             "
                   "  (wall-time: {:.6f}s,  time-step duration: {:.2e}s)\n",
                   100 * s.time.count() / time, s.time.count(), s.dt.count(),
                   s.cycle, wall_time.count(), diff.count());
        write_job = write_job.then(
            [=](::hpx::future<void>) { fub::invoke(feedback, state); });
        start = std::chrono::steady_clock::now();
        return s.time < sub_goal && s.cycle < cycles;
      };
      state = solver.advance(state, state.time + dt, feedback);
    }
    fmt::print("Wait for write jobs to finish...\n");
    write_job.get();
  } catch (std::exception &e) {
    fmt::print(e.what());
  }
}

} // namespace hpx
} // namespace fub

#endif // !DRIVER_HPP
