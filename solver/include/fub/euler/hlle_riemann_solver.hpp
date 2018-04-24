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

#ifndef FUB_EULER_HLLERIEMANNSOLVER_HPP
#define FUB_EULER_HLLERIEMANNSOLVER_HPP

#include "fub/algorithm.hpp"
#include "fub/equation.hpp"
#include "fub/euler/hll_riemann_solver.hpp"
#include "fub/euler/hlle_signal_velocity.hpp"

namespace fub {
namespace euler {

struct hlle_riemann_solver {
  template <typename E> using State = complete_state_t<E>;
  hlle_signal_velocity m_signals{0.5};

  /// @brief Returns the F(q*) where F is the flux defined by the equation, q*
  /// the approximated solution to the Riemann Problem spanned by `left` and
  /// `right` states.
  template <typename Equation>
  add_flux_t<conservative_state_t<Equation>>
  compute_numeric_flux(const State<Equation>& left,
                       const State<Equation>& right,
                       const Equation& equation) const noexcept {
    return hll_riemann_solver::compute_numeric_flux(left, right, m_signals,
                                                    equation);
  }

  /// @brief Returns the F(q*) where F is the flux defined by the equation, q*
  /// the approximated solution to the Riemann Problem spanned by `left` and
  /// `right` states.
  template <typename Abi, typename Equation>
  auto
  compute_numeric_flux(const Abi& abi,
                       const add_simd_t<complete_state_t<Equation>, Abi>& left,
                       const add_simd_t<complete_state_t<Equation>, Abi>& right,
                       const Equation& equation) const noexcept {
    return hll_riemann_solver::compute_numeric_flux(abi, left, right, m_signals,
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
