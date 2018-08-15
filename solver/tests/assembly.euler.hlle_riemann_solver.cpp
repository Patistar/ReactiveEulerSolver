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

constexpr int Rank = 3;

using Equation = fub::euler::perfect_gas<Rank>;
using EquationF = fub::euler::perfect_gas<Rank, float>;
using EquationLD = fub::euler::perfect_gas<Rank, long double>;

void solve_double(Equation equation, fub::const_patch_t<Equation> left,
                  fub::const_patch_t<Equation> middle,
                  fub::const_patch_t<Equation> right,
                  fub::fluxes_t<Equation> flux) {
  fub::euler::hlle_riemann_solver::compute_numeric_fluxes(
      equation, left, middle, right, flux);
}