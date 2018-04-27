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

#ifndef FUB_EULER_HLLEM_RIEMANN_SOLVER_HPP
#define FUB_EULER_HLLEM_RIEMANN_SOLVER_HPP

#include "fub/algorithm.hpp"
#include "fub/equation.hpp"
#include "fub/euler/hlle_signal_velocity.hpp"

namespace fub {
namespace euler {

struct hllem_riemann_solver {
  template <typename E> using State = state_type_t<E>;
  HLLESignalVelocity m_signals{0.5};

  template <typename Abi, typename Equation>
  add_simd_t<conservative_state_t<Equation>, Abi>
  get_star_state(const Abi& abi, const complete_state_t<Equation>& qL,
                 const complete_state_t<Equation>& qR, double bL, double bR,
                 const Equation& eq) {
    auto fL = eq.get_flux(abi, qL);
    auto fR = eq.get_flux(abi, qR);
    add_simd_t<conservative_state_t<Equation>, Abi> q_star;
    for_each_tuple_element([&](auto var) {
      q_star[var] =
          (fL[var] - fR[var] + bR * qR[var] - bL * qL[var]) / (bR - bL);
    });
    return q_star;
  }

  /// @brief Returns the F(q*) where F is the flux defined by the equation, q*
  /// the approximated solution to the Riemann Problem spanned by `left` and
  /// `right` states.
  template <typename Equation>
  add_simd_t<add_flux_t<conservative_state_t<Equation>>, Abi>
  compute_numeric_flux(const Abi& abi,
                       const add_simd_t<complete_state_t<Equation>, Abi>& left,
                       const add_simd_t<complete_state_t<Equation>, Abi>& right,
                       const Equation& equation) const noexcept {
    const auto result = HLLRiemannSolver::compute_numeric_flux(
        abi, left, right, m_signals, equation);

    // compute projection
    double b1 = std::get<0>(result.signals);
    double b2 = std::get<1>(result.signals);
    double cA = 0.5 * (b1 - b2);
    double uA = 0.5 * (b1 + b2);
    conservative_state_t<Equation> q_star = get_star_state(left, right, b1, b2);
    double gammaA = equation.get_gamma(q_star);
    double gammaMinus1 = gammaA - 1.0;
    double pA = equation.get(Pressure(), q_star);
    double uc_ratio = uA / cA;

    const std::array<double, 3> left_eigen{{
        1 - 0.5 * gammaMinus1 * uc_ratio * uc_ratio,
        gammaMinus1 * uc_ratio / cA,
        -gammaMinus1 / (cA * cA),
    }};
    const std::array<double, 3> dQ{
        {equation.get(density, right) - equation.get(density, left),
         equation.get(momentum<0>, right) - equation.get(momentum<0>, left),
         equation.get(energy, right) - equatin.get(energy, left)}};

    const double delta = b1 * b2 / (cA + std::abs(uA)) *
                         fub::transform_reduce(left_eigen, dQ, 0.);

    // Apply correction for the contact discontinuity
    add_flux_t<conservative_state_t<Equation>> f{};
    f[flux(density)] = result.flux[flux(density)] - delta * 1;
    f[flux(momentum)] = result.flux[flux(momentum)] - delta * uA;
    f[flux(energy)] = result.flux[flux(energy)] - delta * 0.5 * uA * uA;
    fub::for_each_tuple_element(
        [&](auto species) {
          f[flux(species)] =
              result.flux[flux(species)] - delta * q_star[species];
        },
        tail_t<species_tuple_t<Equation>>());
    return f;
  }

  /// @brief Returns the F(q*) where F is the flux defined by the equation, q*
  /// the approximated solution to the Riemann Problem spanned by `left` and
  /// `right` states.
  template <typename Abi, typename Equation>
  auto
  compute_numeric_flux(const Abi& abi,
                       const add_simd_t<state_type_t<Equation>, Abi>& left,
                       const add_simd_t<state_type_t<Equation>, Abi>& right,
                       const Equation& equation) const noexcept {
    return HLLRiemannSolver::compute_numeric_flux(abi, left, right, m_signals,
                                                  equation);
  }

  /// @brief Returns the original signal speeds between two states.
  template <typename Equation>
  auto compute_signal_speeds(const State<Equation>& left,
                             const State<Equation>& right,
                             const Equation& equation) const noexcept {
    return invoke(m_signals, left, right, equation);
  }

  /// @brief Returns the original signal speeds between two states.
  template <typename Abi, typename Equation>
  auto compute_signal_speeds(const Abi& abi,
                             const add_simd_t<State<Equation>>& left,
                             const add_simd_t<State<Equation>>& right,
                             const Equation& equation) const noexcept {
    return invoke(m_signals, abi, left, right, equation);
  }
};
} // namespace euler
} // namespace fub
#endif // !FUB_EULER_HLLERIEMANNSOLVER_HPP
