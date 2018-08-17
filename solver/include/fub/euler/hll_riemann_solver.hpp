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

#ifndef FUB_EULER_HLLRIEMANNSOLVER_HPP
#define FUB_EULER_HLLRIEMANNSOLVER_HPP

#include "fub/equation.hpp"
#include "fub/euler/variables.hpp"

namespace fub {
inline namespace v1 {
namespace euler {

struct hll_riemann_solver {
  template <typename SignalFunction, typename Equation, typename StateRefL,
            typename StateRefR, typename FluxRef>
  static void do_compute_numeric_flux(SignalFunction signal_fn,
                                      Equation equation, const StateRefL& qL,
                                      const StateRefR& qR, const FluxRef& fL, 
                                      const FluxRef& fR,
                                      const FluxRef& flux) {
    using float_v = typename FluxRef::value_type;
    auto signals = signal_fn(equation, qL, qR);
    float_v sL = signals.left;
    float_v sR = signals.right;
    equation.flux(qL, fL);
    equation.flux(qR, fR);
    for_each(flux.get_variable_list(), [&](auto var) {
      float_v f =
          (sR * fL[var] - sL * fR[var] + sL * sR * (qR[var] - qL[var])) /
          (sR - sL);
      where(sL >= 0, f) = fL[var];
      where(sR <= 0, f) = fR[var];
      flux[var] = f;
    });
  }

  template <typename SignalFunction, typename Equation, typename StateRefL,
            typename StateRefR, typename FluxRef>
  static void do_compute_numeric_flux(SignalFunction signal_fn,
                                      Equation equation, const StateRefL& qL,
                                      const StateRefR& qR,
                                      const FluxRef& flux) {
    flux_state_t<Equation> fL(flux.get_variable_list());
    flux_state_t<Equation> fR(flux.get_variable_list());
    do_compute_numeric_flux(signal_fn, equation, qL, qR, fL(0), fR(0), flux);
  }

  template <typename SignalFunction, typename Equation>
  static void compute_numeric_fluxes(SignalFunction signal_fn,
                                     Equation equation,
                                     const_patch_t<Equation> left,
                                     const_patch_t<Equation> mid,
                                     const_patch_t<Equation> right,
                                     fluxes_t<Equation> fluxes) {
    for_each_index(mid.get_mapping(), [&](auto... is) {
      constexpr int Rank = Equation::rank();
      std::array<std::ptrdiff_t, Rank> iL{{is...}};
      std::array<std::ptrdiff_t, Rank> iR = shift(iL, 0, 1);
      if (iL[0] == 0) {
        const int last_extent_left = left.get_extents().extent(0) - 1;
        std::array<std::ptrdiff_t, Rank> iLL = replace(iL, 0, last_extent_left);
        do_compute_numeric_flux(signal_fn, equation, left(iLL), mid(iL),
                                fluxes(iL));
      }
      if (iR[0] < mid.get_extents().extent(0)) {
        do_compute_numeric_flux(signal_fn, equation, mid(iL), mid(iR),
                                fluxes(iR));
      } else {
        std::array<std::ptrdiff_t, Rank> iR_ = replace(iR, 0, 0);
        do_compute_numeric_flux(signal_fn, equation, mid(iL), right(iR_),
                                fluxes(iR));
      }
    });
  }
};

} // namespace euler
}
} // namespace fub

#endif
