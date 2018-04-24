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

#include "fub/flux_method.hpp"
#include "fub/patch_view.hpp"
#include "fub/simd.hpp"
#include "fub/utility.hpp"

#include <chrono>
#include <iterator>

namespace fub {

template <typename RiemannSolver> struct godunov_method {
  RiemannSolver m_riemann_solver;

  template <typename Abi, typename Equation, typename Coordinates, typename L,
            typename R, typename F>
  std::enable_if_t<is_simd_abi_v<Abi>>
  compute_numeric_fluxes(const Abi &abi, const F &flux, const L &left,
                         const R &right, std::chrono::duration<double>,
                         const Coordinates &, const Equation &equation) const {
    static_assert(
        conjunction<is_view<L>, is_view<R>, is_view<F>>::value,
        "compute_numeric_fluxes was not invoked with a patch_view type. "
        "Consider using the make_view function.");
    for_each_simd(abi,
                  [&](auto abi, auto &&f, const auto &qL, const auto &qR) {
                    f = m_riemann_solver.compute_numeric_flux(abi, qL, qR,
                                                              equation);
                  },
                  flux, left, right);
  }

  template <typename Abi, typename Equation, typename Coordinates, typename L,
            typename M, typename R, typename F>
  std::enable_if_t<is_simd_abi_v<Abi>> compute_numeric_fluxes(
      const Abi &abi, const F &flux, const L &left, const M &middle,
      const R &right, std::chrono::duration<double> dt,
      const Coordinates &coords, const Equation &equation) const {
    static_assert(
        conjunction<is_view<L>, is_view<M>, is_view<R>, is_view<F>>::value,
        "compute_numeric_fluxes was not invoked with a patch_view type. "
        "Consider using the make_view function.");
    for_each_row(
        [&](auto left, auto middle, auto right, auto flux) {
          const auto row = join(rtake<1>(left), middle, take<1>(right));
          auto qL = rdrop<1>(make_view(row));
          auto qR = drop<1>(make_view(row));
          compute_numeric_fluxes(abi, flux, qL, qR, dt, coords, equation);
        },
        left, middle, right, flux);
  }

  template <typename Equation, typename Coordinates, typename L, typename M,
            typename R, typename F>
  std::enable_if_t<
      conjunction<is_view<F>, is_view<L>, is_view<M>, is_view<R>>::value>
  compute_numeric_fluxes(const F &flux, const L &left, const M &middle,
                         const R &right, std::chrono::duration<double> dt,
                         const Coordinates &coordinates,
                         const Equation &equation) const {
    return compute_numeric_fluxes(simd_abi::native<double>(), flux, left,
                                  middle, right, dt, coordinates, equation);
  }

  template <typename Equation, typename L, typename M, typename R,
            typename Coordinates>
  std::enable_if_t<is_view<L>::value && is_view<R>::value && is_view<M>::value,
                   std::chrono::duration<double>>
  get_stable_time_step(const Equation &equation, const L &left, const M &middle,
                       const R &right, const Coordinates &coordinates) const {
    double max_S{0};
    for_each_row(
        [&](auto left, auto middle, auto right) {
          auto max = [](auto current, const auto &tuple) {
            auto all_signals = std::tuple_cat(std::make_tuple(current), tuple);
            return fub::apply(
                [](auto... signals) {
                  auto abs = [](auto x) { return x < 0 ? -x : +x; };
                  return std::max({abs(signals)...});
                },
                all_signals);
          };
          max_S = max(max_S, m_riemann_solver.compute_signal_speeds(
                                 left.last(), middle.first(), equation));
          max_S = max(max_S, m_riemann_solver.compute_signal_speeds(
                                 middle.last(), right.first(), equation));
          auto it = middle.begin();
          const auto last = middle.end();
          while (it != last && std::next(it) != last) {
            max_S = max(max_S, m_riemann_solver.compute_signal_speeds(
                                   *it, *std::next(it), equation));
            ++it;
          }
        },
        left, middle, right);
    assert(max_S > 0);
    double dx = coordinates.dx()[0];
    return std::chrono::duration<double>(0.5 * dx / max_S);
  }
};

} // namespace fub

#endif
