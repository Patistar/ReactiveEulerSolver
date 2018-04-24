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
#include "fub/patch_view.hpp"

#include <cassert>
#include <chrono>
#include <functional>

namespace fub {
namespace euler {

class hll_riemann_solver {
  static double apply_formula(double uL, double uR, double fL, double fR,
                              double sL, double sR) noexcept;

  static scalar<double>
  apply_formula(simd_abi::scalar, const scalar<double>& uL,
                const scalar<double>& uR, const scalar<double>& fL,
                const scalar<double>& fR, const scalar<double>& sL,
                const scalar<double>& sR) noexcept;

#if defined(FUB_SIMD_HAS_SSE)
  static sse<double> apply_formula(simd_abi::sse, const sse<double>& uL,
                                   const sse<double>& uR, const sse<double>& fL,
                                   const sse<double>& fR, const sse<double>& sL,
                                   const sse<double>& sR) noexcept;
#endif

#if defined(FUB_SIMD_HAS_AVX)
  static avx<double> apply_formula(simd_abi::avx, const avx<double>& uL,
                                   const avx<double>& uR, const avx<double>& fL,
                                   const avx<double>& fR, const avx<double>& sL,
                                   const avx<double>& sR) noexcept;
#endif

  /// @brief Returns the HLL flux for all variables of the specified `Flux`
  /// type.
  ///
  /// @note This function assumes, that equation fluxes for left and right
  /// states are already given.
  template <typename Signal, typename Equation>
  static flux_type_t<Equation> do_compute_numeric_flux(
      Signal signal_fn, const complete_state_t<Equation>& qL,
      const flux_type_t<Equation>& fL, const complete_state_t<Equation>& qR,
      const flux_type_t<Equation>& fR, const Equation& equation) noexcept {
    const auto signals = invoke(signal_fn, qL, qR, equation);
    const auto compute = [&](auto quantity) {
      return apply_formula(qL[quantity], qR[quantity], fL[flux(quantity)],
                           fR[flux(quantity)], std::get<0>(signals),
                           std::get<1>(signals));
    };
    return fub::apply(
        [&](auto... vars) { return flux_type_t<Equation>{compute(vars)...}; },
        get_variables(conservative_state_t<Equation>()));
  }

  template <typename Equation, typename Abi>
  using simd_state_t = add_simd_t<complete_state_t<Equation>, Abi>;

  template <typename Equation, typename Abi>
  using simd_flux_t =
      add_flux_t<add_simd_t<conservative_state_t<Equation>, Abi>>;

  /// @brief Returns the HLL flux for all variables of the specified `Flux`
  /// type.
  ///
  /// @note This function assumes, that equation fluxes for left and right
  /// states are already given.
  template <typename Abi, typename Signal, typename Equation>
  static simd_flux_t<Equation, Abi> do_compute_numeric_flux_simd(
      const Abi& abi, Signal signal_fn, const simd_state_t<Equation, Abi>& qL,
      const simd_flux_t<Equation, Abi>& fL,
      const simd_state_t<Equation, Abi>& qR,
      const simd_flux_t<Equation, Abi>& fR, const Equation& equation) noexcept {
    const auto signals = invoke(signal_fn, abi, qL, qR, equation);
    const auto compute = [&](auto quantity) {
      return apply_formula(abi, qL[quantity], qR[quantity], fL[flux(quantity)],
                           fR[flux(quantity)], std::get<0>(signals),
                           std::get<1>(signals));
    };
    return fub::apply(
        [&](auto... vars) {
          return simd_flux_t<Equation, Abi>{compute(vars)...};
        },
        get_variables(conservative_state_t<Equation>()));
  }

public:
  /// @brief Returns the HLL flux between to states `left` and `right`.
  ///
  /// @detail This function uses the specified equation to compute state fluxes
  template <typename Signal, typename Equation>
  static flux_type_t<Equation>
  compute_numeric_flux(const complete_state_t<Equation>& left,
                       const complete_state_t<Equation>& right,
                       Signal signal_fn, const Equation& equation) noexcept {
    const flux_type_t<Equation> flux_left = equation.get_flux(left);
    const flux_type_t<Equation> flux_right = equation.get_flux(right);
    return hll_riemann_solver::do_compute_numeric_flux(
        signal_fn, left, flux_left, right, flux_right, equation);
  }

  template <typename Abi, typename Signal, typename Equation>
  static simd_flux_t<Equation, Abi>
  compute_numeric_flux(const Abi& abi,
                       const add_simd_t<complete_state_t<Equation>, Abi>& left,
                       const add_simd_t<complete_state_t<Equation>, Abi>& right,
                       Signal signal_fn, const Equation& equation) noexcept {
    const simd_flux_t<Equation, Abi> flux_left = equation.get_flux(abi, left);
    const simd_flux_t<Equation, Abi> flux_right = equation.get_flux(abi, right);
    return hll_riemann_solver::do_compute_numeric_flux_simd(
        abi, signal_fn, left, flux_left, right, flux_right, equation);
  }
};

} // namespace euler
} // namespace fub

#endif
