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

#ifndef FUB_SOLVER_GODUNOV_METHOD_HPP
#define FUB_SOLVER_GODUNOV_METHOD_HPP

#include "fub/equation.hpp"

#include <chrono>

namespace fub {

template <typename RiemannSolver> struct godunov_method {
  // static_assert(is_riemann_solver<RiemannSolver>(),
  //               "RiemannSolver does not fulfill the RiemannSolver concept");

  RiemannSolver m_riemann_solver;

  template <typename T> using duration_t = std::chrono::duration<T>;

  template <typename Equation, typename T>
  void compute_numeric_fluxes(Equation equation, duration_t<T> /* dt */,
                              T /* dx */, const_patch_t<Equation> left,
                              const_patch_t<Equation> mid,
                              const_patch_t<Equation> right,
                              fluxes_t<Equation> fluxes) const {
    m_riemann_solver.compute_numeric_fluxes(equation, left, mid, right, fluxes);
  }

  template <typename Equation, typename T>
  duration_t<T> estimate_time_step_size(Equation equation, T dx,
                                        const_patch_t<Equation> left,
                                        const_patch_t<Equation> mid,
                                        const_patch_t<Equation> right) const {
    T max_wave_speed = 0;
    for_each_face(left, mid, right, [&](auto qL, auto qR) {
      auto signals = m_riemann_solver.compute_signals(equation, qL, qR);
      max_wave_speed = std::max({max_wave_speed, -signals.left, signals.right});
    });
    assert(max_wave_speed > 0.);
    assert(dx > 0.);
    return dx / max_wave_speed;
  }
};

} // namespace fub

#endif
