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

#ifndef FUB_ODE_SOLVER_ROSENBROCK_HPP
#define FUB_ODE_SOLVER_ROSENBROCK_HPP

#include <boost/numeric/odeint.hpp>

namespace fub {
namespace ode_solver {

struct rosenbrock4 {
  using vector_type = boost::numeric::ublas::vector<double>;
  using matrix_type = boost::numeric::ublas::matrix<double>;

  template <typename System, index Size>
  static void integrate(const System& system, span<double, Size> y0,
                        std::chrono::duration<double> dt) {
    auto system_blas = [system](const vector_type& x, vector_type& dxdt,
                                double t) {
      span<const double, Size> x_(x.data());
      span<double, Size> dxdt_(dxdt.data());
      system(dxdt_, x_, t);
    };
    auto system_jacobi = [system_blas](const vector_type& y0, matrix_type& Dy,
                                       double x, vector_type& dfdt) {
      static constexpr double eps = 1.0e-16;
      std::fill(dfdt.begin(), dfdt.end(), 0.0);
      vector_type dydx(Size);
      vector_type dydx_0(Size);
      vector_type y(y0);
      system_blas(y0, dydx_0, x);
      for (index i = 0; i < Size; ++i) {
        const double delta =
            std::sqrt(eps * std::max(1.0e-5, std::fabs(y0[i])));
        y[i] += delta;
        system_blas(y, dydx, x);
        for (index j = 0; j < Size; ++j) {
          Dy(j, i) = (dydx[j] - dydx_0[j]) / delta;
        }
        y[i] = y0[i];
      } 
    };
    vector_type start_state(Size);
    std::copy(y0.begin(), y0.end(), start_state.begin());
    boost::numeric::odeint::integrate_adaptive(
        boost::numeric::odeint::rosenbrock4<double>{},
        std::make_pair(system_blas, system_jacobi), start_state, 0., dt.count(),
        dt.count());
    std::copy(start_state.begin(), start_state.end(), y0.begin());
  }

  template <typename Archive> void serialize(Archive&, unsigned) {}
};

} // namespace ode_solver
} // namespace fub

#endif // FUB_ODE_SOLVER_ROSENBROCK_HPP