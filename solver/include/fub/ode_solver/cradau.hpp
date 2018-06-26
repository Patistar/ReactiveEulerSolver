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

#ifndef FUB_ODE_SOLVER_CRADAU_HPP
#define FUB_ODE_SOLVER_CRADAU_HPP

#include "fub/functional.hpp"
#include "fub/optional.hpp"
#include "fub/span.hpp"

#include <chrono>
#include <limits>
#include <vector>

namespace fub {
namespace ode_solver {

struct radau5 {
  struct integration_info {
    index newton_iteration_count;
  };

  using ode_system_t =
      function_ref<void(double, span<const double>, span<double>)>;

  using feedback_t =
      function_ref<void(double, span<const double>, const integration_info&)>;

  template <typename System, index Size>
  void integrate(System&& rhs, span<double, Size> y0,
                 std::chrono::duration<double> dt) const {
    auto radau_rhs = [&](double t, span<const double> y_,
                         span<double> dydt_) noexcept {
      span<double, Size> dydt = take<Size>(dydt_);
      span<const double, Size> y = take<Size>(y_);
      fub::invoke(rhs, dydt, y, t);
    };
    std::vector<char> buffer(sizeof(double) * (5 * Size * Size + 12 * Size) +
                             sizeof(int) * 2 * Size + 64 * 1000);
    integrate(radau_rhs, {0, dt.count()}, y0,
              span<char>(buffer.data(), buffer.size()));
  }

  void integrate(ode_system_t ode_system, std::array<double, 2> x,
                 span<double> y_0, span<char> memory,
                 optional<feedback_t> feedback = nullopt) const;

  double absolute_tolerance{1.0E-4};
  double relative_tolerance{1.0E-4};
  double initial_step_size{1.0E-6};
  index max_newton_iteration_count{10};
  int max_restart_count{20};
};

} // namespace ode_solver
} // namespace fub

#endif