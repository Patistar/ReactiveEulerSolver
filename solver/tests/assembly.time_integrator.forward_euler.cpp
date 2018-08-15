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

#include "fub/euler/hlle_riemann_solver.hpp"
#include "fub/euler/perfect_gas.hpp"
#include "fub/godunov_method.hpp"
#include "fub/time_integrator/forward_euler.hpp"

constexpr int Rank = 3;

using Equation = fub::euler::perfect_gas<Rank>;
using FluxMethod = fub::godunov_method<fub::euler::hlle_riemann_solver>;

void integrate_double(const Equation& equation, const FluxMethod& method,
                      std::chrono::duration<double> dt, double dx,
                      fub::const_patch_t<Equation> left,
                      fub::const_patch_t<Equation> mid,
                      fub::const_patch_t<Equation> right,
                      fub::patch_t<Equation> next,
                      fub::fluxes_t<Equation> fluxes) {
  fub::time_integrator::forward_euler integrator{};
  integrator.integrate(equation, method, dt, dx, left, mid, right, next,
                       fluxes);
}