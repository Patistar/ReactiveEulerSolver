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

#include "fub/ode_solver/cradau.hpp"
#include "fub/ode_solver/radau.hpp"
#include "fub/span.hpp"

#include <fmt/format.h>
#include <fmt/ostream.h>

int main() {
  fub::index counter{0};
  auto ode_system = [&](double /* x */, fub::span<const double> y,
                        fub::span<double> dydx) {
    static constexpr double eps = 1.0e-6;
    dydx[0] = y[1];
    dydx[1] = ((1 - y[0] * y[0]) * y[1] - y[0]) / eps;
    ++counter;
  };
  auto print = [&](double x, fub::span<const double> y,
                   const fub::ode_solver::radau5::integration_info&) {
    fmt::print("{} {} {} {}\n", counter, x, y[0], y[1]);
  };
  std::vector<char> buffer(1000);
  std::array<double, 2> y{{2., 0.}};
  fub::ode_solver::radau5 radau{};
  radau.initial_step_size = 1E-6;
  radau.max_newton_iteration_count = 10;
  radau.max_restart_count = 20;
  radau.absolute_tolerance = 1e-5;
  radau.relative_tolerance = 1e-4;
  radau.integrate(ode_system, {0, 11.}, fub::make_span(y), buffer);
}