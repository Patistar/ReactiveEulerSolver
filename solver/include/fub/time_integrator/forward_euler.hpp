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

#include "fub/algorithm.hpp"
#include "fub/equation.hpp"
#include "fub/face.hpp"
#include "fub/patch.hpp"
#include "fub/patch_view.hpp"
#include "fub/simd.hpp"
#include "fub/utility.hpp"

#include <chrono>

namespace fub {
namespace time_integrator {

struct forward_euler {
private:
  static double do_integrate_impl(double lambda, double u, double fL,
                                  double fR) noexcept;

  static scalar<double> do_integrate_impl(simd_abi::scalar,
                                          scalar<double> lambda,
                                          scalar<double> u, scalar<double> fL,
                                          scalar<double> fR) noexcept;

#if defined(FUB_SIMD_HAS_SSE)
  static simd<double> do_integrate_impl(simd_abi::native<double>,
                                        simd<double> lambda, simd<double> u,
                                        simd<double> fL,
                                        simd<double> fR) noexcept;
#endif

  template <typename StatesIn, typename Flux, typename StatesOut,
            typename Equation>
  static void do_integrate(double lambda, StatesIn u, Flux f, StatesOut next,
                           const Equation &equation) noexcept {
    fub::for_each_row(
        [&](auto u, auto f, auto next) {
          auto fL = rdrop<1>(f);
          auto fR = drop<1>(f);
          fub::for_each_simd(
              [&](auto abi, auto &&next, auto u, const auto &fL,
                  const auto &fR) {
                fub::for_each_tuple_element(
                    [&](auto q) {
                      const auto flux_q = flux(q);
                      u[q] = do_integrate_impl(abi, lambda, u[q], fL[flux_q],
                                               fR[flux_q]);
                    },
                    conservative_variables_t<Equation>());
                next = equation.get_complete_state(abi, u);
              },
              next, u, fL, fR);
        },
        u, f, next);
  }

  template <typename Out, typename L, typename M, typename R,
            typename Coordinates, typename FluxMethod, typename Equation>
  static void integrate_impl(const Out &out, const L &left, const M &middle,
                             const R &right, std::chrono::duration<double> dt,
                             const Coordinates &coordinates,
                             const FluxMethod &flux_method,
                             const Equation &equation,
                             std::integral_constant<axis, axis::x>) {

    auto fluxes =
        forward_euler::make_flux_patch<axis::x, Equation>(middle.extents());
    flux_method.compute_numeric_fluxes(make_view(fluxes), left, middle, right,
                                       dt, coordinates, equation);
    double lambda = dt.count() / coordinates.dx()[0];
    do_integrate(lambda, middle, make_view(fub::as_const(fluxes)), out,
                 equation);
  }

  template <axis Axis, typename Out, typename L, typename M, typename R,
            typename Coordinates, typename FluxMethod, typename Equation>
  static void
  integrate_impl(const Out &out, const L &left, const M &middle, const R &right,
                 std::chrono::duration<double> dt,
                 const Coordinates &coordinates, const FluxMethod &flux_method,
                 const Equation &equation, std::integral_constant<axis, Axis>) {
    static constexpr int Dim = as_int(Axis);
    const auto sliced_left = slice_left<Dim, 2>(left);
    const auto permutated_left = permutate<Dim, 0>(make_view(sliced_left));
    const auto sliced_right = slice_right<Dim, 2>(right);
    const auto permutated_right = permutate<Dim, 0>(make_view(sliced_right));
    auto permutated_mid = permutate<Dim, 0>(middle);
    auto fluxes = forward_euler::make_flux_patch<axis::x, Equation>(
        permutated_mid.extents());
    flux_method.compute_numeric_fluxes(
        make_view(fluxes), make_view(permutated_left),
        make_view(fub::as_const(permutated_mid)), make_view(permutated_right),
        dt, coordinates, equation);
    double lambda = dt.count() / coordinates.dx()[as_int(Axis)];
    do_integrate(lambda, make_view(as_const(permutated_mid)),
                 make_view(as_const(fluxes)), make_view(permutated_mid),
                 equation);
    permutate<as_int(Axis), 0>(out, make_view(as_const(permutated_mid)));
  }

public:
  template <axis Axis, typename Out, typename L, typename M, typename R,
            typename Coordinates, typename FluxMethod, typename Equation>
  static void
  integrate(const Out &out, const L &left, const M &middle, const R &right,
            std::chrono::duration<double> dt, const Coordinates &coordinates,
            const FluxMethod &flux_method, const Equation &equation) {
    integrate_impl(out, left, middle, right, dt, coordinates, flux_method,
                   equation, std::integral_constant<axis, Axis>());
  }

  template <axis Axis, typename Equation, typename Extents>
  static auto make_flux_patch(const Extents &extents) {
    static constexpr int Dim = as_int(Axis);
    auto grown = grow(extents, int_c<Dim>);
    return make_patch(as_tuple_t<flux_type_t<Equation>>(), grown);
  }
};

} // namespace time_integrator
} // namespace fub

#endif
