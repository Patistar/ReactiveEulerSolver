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

#include "fub/equation.hpp"
#include "fub/euler/hll_riemann_solver.hpp"
#include "fub/euler/hlle_signal_velocity.hpp"

namespace fub {
inline namespace v1 {
namespace euler {

/// \ingroup RiemannSolver
/// This class computes the numeric flux via is the approximate HLL Riemann
/// Solver with Einfeldt's signal velocities.
///
/// \note This does not apply Einfeldts correction for the moving contact
/// discontinuity.
struct hlle_riemann_solver {
  /// Computes F(q*) where F is the flux defined by the equation, q*
  /// the approximated solution to the Riemann Problem spanned by `left` and
  /// `right` states.
  ///
  /// \param[in] equation  The base equation which knows of all formulas.
  /// \param[in] left  The left face neighbor patch of middle.
  /// \param[in] middle  This is the patch which the fluxes are computed for.
  /// \param[in] right  The right face neighbor patch of middle.
  /// \param[out] fluxes  A face patch where the results will be stored at.
  template <typename Equation>
  static void compute_numeric_fluxes(Equation equation,
                                     const_patch_t<Equation> left,
                                     const_patch_t<Equation> middle,
                                     const_patch_t<Equation> right,
                                     fluxes_t<Equation> flux) noexcept {
    return hll_riemann_solver::compute_numeric_fluxes(
        hlle_signal_velocity{}, equation, left, middle, right, flux);
  }
};

} // namespace euler
} // namespace v1
} // namespace fub

#endif // !FUB_EULER_HLLERIEMANNSOLVER_HPP
