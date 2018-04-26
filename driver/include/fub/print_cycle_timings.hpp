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

#ifndef FUB_DRIVER_PRINT_CYCLE_TIMINGS_HPP
#define FUB_DRIVER_PRINT_CYCLE_TIMINGS_HPP

#include <chrono>
#include <fmt/format.h>

namespace fub {

class print_cycle_timings {
  std::optional<std::chrono::duration<double>> final_time{};
  std::chrono::steady_clock::time_point wall{std::chrono::steady_clock::now()};
  std::chrono::steady_clock::time_point cycle{std::chrono::steady_clock::now()};

public:
  print_cycle_timings() = default;

  explicit print_cycle_timings(std::chrono::duration<double> t_end)
      : final_time{t_end} {}

  template <typename State> bool operator()(const State& s) {
    auto end = std::chrono::steady_clock::now();
    auto wall_time =
        std::chrono::duration_cast<std::chrono::duration<double>>(end - wall);
    auto cycle_time =
        std::chrono::duration_cast<std::chrono::duration<double>>(end - cycle);
    if (final_time) {
      fmt::print("[{:.2f}%] T = {:.10f}s, dt = {:.4e}s, cycle = {:10}   "
                 "             "
                 "  (wall-time: {:.6f}s,  time-step duration: {:.2e}s)\n",
                 100 * (s.time / *final_time).get(), s.time.count(),
                 s.dt.count(), s.cycle, wall_time.count(), cycle_time.count());
    } else {
      fmt::print("T = {:.10f}s, dt = {:.4e}s, cycle = {:10}   "
                 "             "
                 "  (wall-time: {:.6f}s,  time-step duration: {:.2e}s)\n",
                 s.time.count(), s.dt.count(), s.cycle, wall_time.count(),
                 cycle_time.count());
    }
    cycle = std::chrono::steady_clock::now();
    return true;
  }
};

} // namespace fub

#endif // !PRINT_CYCLE_TIMINGS_HPP
