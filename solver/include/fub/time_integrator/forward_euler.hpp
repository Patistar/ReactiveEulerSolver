// Copyright (c) 2017 Maikel Nadolski
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

#ifndef FUB_FORWARD_EULER_HPP
#define FUB_FORWARD_EULER_HPP

#include "fub/equation.hpp"

#include <chrono>

namespace fub {
inline namespace v1 {
namespace time_integrator {

struct forward_euler {
  template <typename T> using duration_t = std::chrono::duration<T>;

  template <typename Equation, typename T, typename FluxMethod>
  static void integrate(const Equation& equation, const FluxMethod& method,
                        duration_t<T> dt, T dx, const_patch_t<Equation> left,
                        const_patch_t<Equation> mid,
                        const_patch_t<Equation> right, 
                        patch_t<Equation> next,
                        fluxes_t<Equation> fluxes) {
    static constexpr int Rank = Equation::rank();
    method.compute_numeric_fluxes(equation, dt, dx, left, mid, right, fluxes);
    const double lambda = dt.count() / dx;
    for_each_index(next.get_mapping(), [&](auto... indices) {
      std::array<index, Rank> left{{indices...}};
      std::array<index, Rank> right = shift(left, 0, 1);
      for_each(fluxes.get_variable_list(), [&](auto var) {
        next[var](indices...) =
            mid[var](indices...) +
            lambda * (fluxes[var](left) - fluxes[var](right));
      });
      equation.reconstruct(next(indices...), next(indices...));
    });
  }
};

} // namespace time_integrator
} // namespace v1
} // namespace fub

#endif
