// Copyright (c) 2018 Maikel Nadolski
//
// Permission is hereby granted, free of chArge, to any person obtaining a copy
// of this softwAre and associated documentation files (the "Software"), to deal
// in the SoftwAre without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the SoftwAre, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the SoftwAre.
//
// THE SOFTWArE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WArRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PArTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ArISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWArE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWArE.

#include "fub/euler/mechanism/burke_2012.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace fub {
namespace euler {
namespace mechanism {
namespace burke2012 {

template <typename T> using span = Burke2012::species_span<T>;

template <typename Species>
static constexpr int index_v =
    variable_find_index<Species, Burke2012::species_tuple>::value;

static double W[std::tuple_size<Burke2012::species_tuple>::value];

auto make_molar_masses() noexcept {
  W[index_v<N2>] = 2.80200000e+01;
  W[index_v<Ar>] = 3.99480000e+01;
  W[index_v<H>] = 1.00800000e+00;
  W[index_v<O2>] = 3.20000000e+01;
  W[index_v<O>] = 1.60000000e+01;
  W[index_v<OH>] = 1.70080000e+01;
  W[index_v<H2>] = 2.01600000e+00;
  W[index_v<H2O>] = 1.80160000e+01;
  W[index_v<HE>] = 4.00000000e+00;
  W[index_v<HO2>] = 3.30080000e+01;
  W[index_v<H2O2>] = 3.40160000e+01;
  return W;
}

span<const double> Burke2012::get_molar_masses() noexcept {
  static auto molar_masses = make_molar_masses();
  return span<const double>{&molar_masses[0],
                            std::tuple_size<Burke2012::species_tuple>::value};
}

namespace {
template <typename Arithmetic>
std::enable_if_t<std::is_arithmetic<Arithmetic>::value>
get_specific_heat_capacities_at_constant_pressure_impl(span<Arithmetic> cp,
                                                       Arithmetic T) noexcept {
  T = std::max(T, Arithmetic(300));
  if (T > 1000.0) {
    cp[index_v<N2>] = 2.96733125e+02 *
                      (2.92664000e+00 +
                       T * (1.48797700e-03 +
                            T * (-5.68476100e-07 +
                                 T * (1.00970400e-10 + T * -6.75335100e-15))));
    cp[index_v<Ar>] = 2.08132126e+02 *
                      (2.50000000e+00 +
                       T * (0.00000000e+00 +
                            T * (0.00000000e+00 +
                                 T * (0.00000000e+00 + T * 0.00000000e+00))));
    cp[index_v<H>] = 8.24847438e+03 *
                     (2.50000000e+00 +
                      T * (0.00000000e+00 +
                           T * (0.00000000e+00 +
                                T * (0.00000000e+00 + T * 0.00000000e+00))));
    cp[index_v<O2>] = 2.59826943e+02 *
                      (3.69757800e+00 +
                       T * (6.13519700e-04 +
                            T * (-1.25884200e-07 +
                                 T * (1.77528100e-11 + T * -1.13643500e-15))));
    cp[index_v<O>] = 5.19653886e+02 *
                     (2.54206000e+00 +
                      T * (-2.75506200e-05 +
                           T * (-3.10280300e-09 +
                                T * (4.55106700e-12 + T * -4.36805200e-16))));
    cp[index_v<OH>] = 4.88855960e+02 *
                      (2.86472886e+00 +
                       T * (1.05650448e-03 +
                            T * (-2.59082758e-07 +
                                 T * (3.05218674e-11 + T * -1.33195876e-15))));
    cp[index_v<H2>] = 4.12423719e+03 *
                      (2.99142300e+00 +
                       T * (7.00064400e-04 +
                            T * (-5.63382900e-08 +
                                 T * (-9.23157800e-12 + T * 1.58275200e-15))));
    cp[index_v<H2O>] = 4.61504339e+02 *
                       (2.67214600e+00 +
                        T * (3.05629300e-03 +
                             T * (-8.73026000e-07 +
                                  T * (1.20099600e-10 + T * -6.39161800e-15))));
    cp[index_v<HE>] = 2.07861554e+03 *
                      (2.50000000e+00 +
                       T * (0.00000000e+00 +
                            T * (0.00000000e+00 +
                                 T * (0.00000000e+00 + T * 0.00000000e+00))));

    cp[index_v<HO2>] = 2.51892334e+02 *
                       (4.01721090e+00 +
                        T * (2.23982013e-03 +
                             T * (-6.33658150e-07 +
                                  T * (1.14246370e-10 + T * -1.07908535e-14))));
    cp[index_v<H2O2>] =
        2.44427980e+02 *
        (4.57316700e+00 +
         T * (4.33613600e-03 +
              T * (-1.47468900e-06 +
                   T * (2.34890400e-10 + T * -1.43165400e-14))));
  } else {
    cp[index_v<N2>] = 2.96733125e+02 *
                      (3.29867700e+00 +
                       T * (1.40824000e-03 +
                            T * (-3.96322200e-06 +
                                 T * (5.64151500e-09 + T * -2.44485500e-12))));
    cp[index_v<Ar>] = 2.08132126e+02 *
                      (2.50000000e+00 +
                       T * (0.00000000e+00 +
                            T * (0.00000000e+00 +
                                 T * (0.00000000e+00 + T * 0.00000000e+00))));
    cp[index_v<H>] = 8.24847438e+03 *
                     (2.50000000e+00 +
                      T * (0.00000000e+00 +
                           T * (0.00000000e+00 +
                                T * (0.00000000e+00 + T * 0.00000000e+00))));
    cp[index_v<O2>] = 2.59826943e+02 *
                      (3.21293600e+00 +
                       T * (1.12748600e-03 +
                            T * (-5.75615000e-07 +
                                 T * (1.31387700e-09 + T * -8.76855400e-13))));
    cp[index_v<O>] = 5.19653886e+02 *
                     (2.94642900e+00 +
                      T * (-1.63816600e-03 +
                           T * (2.42103200e-06 +
                                T * (-1.60284300e-09 + T * 3.89069600e-13))));
    cp[index_v<OH>] = 4.88855960e+02 *
                      (4.12530561e+00 +
                       T * (-3.22544939e-03 +
                            T * (6.52764691e-06 +
                                 T * (-5.79853643e-09 + T * 2.06237379e-12))));
    cp[index_v<H2>] = 4.12423719e+03 *
                      (3.29812400e+00 +
                       T * (8.24944200e-04 +
                            T * (-8.14301500e-07 +
                                 T * (-9.47543400e-11 + T * 4.13487200e-13))));
    cp[index_v<H2O>] = 4.61504339e+02 *
                       (3.38684200e+00 +
                        T * (3.47498200e-03 +
                             T * (-6.35469600e-06 +
                                  T * (6.96858100e-09 + T * -2.50658800e-12))));
    cp[index_v<HE>] = 2.07861554e+03 *
                      (2.50000000e+00 +
                       T * (0.00000000e+00 +
                            T * (0.00000000e+00 +
                                 T * (0.00000000e+00 + T * 0.00000000e+00))));
    cp[index_v<HO2>] = 2.51892334e+02 *
                       (4.30179801e+00 +
                        T * (-4.74912051e-03 +
                             T * (2.11582891e-05 +
                                  T * (-2.42763894e-08 + T * 9.29225124e-12))));
    cp[index_v<H2O2>] =
        2.44427980e+02 *
        (3.38875400e+00 +
         T * (6.56922600e-03 +
              T * (-1.48501300e-07 +
                   T * (-4.62580600e-09 + T * 2.47151500e-12))));
  }
}

template <typename Simd>
std::enable_if_t<is_simd<Simd>::value>
get_specific_heat_capacities_at_constant_pressure_impl(span<Simd> cp,
                                                       Simd T) noexcept {
  T = fub::max(T, Simd(300));
  auto big_T = T > 1000.0;
  if (any_of(big_T)) {
    where(big_T, cp[index_v<N2>]) =
        2.96733125e+02 *
        (2.92664000e+00 +
         T * (1.48797700e-03 +
              T * (-5.68476100e-07 +
                   T * (1.00970400e-10 + T * -6.75335100e-15))));
    where(big_T, cp[index_v<Ar>]) =
        2.08132126e+02 *
        (2.50000000e+00 +
         T * (0.00000000e+00 +
              T * (0.00000000e+00 +
                   T * (0.00000000e+00 + T * 0.00000000e+00))));
    where(big_T, cp[index_v<H>]) =
        8.24847438e+03 *
        (2.50000000e+00 +
         T * (0.00000000e+00 +
              T * (0.00000000e+00 +
                   T * (0.00000000e+00 + T * 0.00000000e+00))));
    where(big_T, cp[index_v<O2>]) =
        2.59826943e+02 *
        (3.69757800e+00 +
         T * (6.13519700e-04 +
              T * (-1.25884200e-07 +
                   T * (1.77528100e-11 + T * -1.13643500e-15))));
    where(big_T, cp[index_v<O>]) =
        5.19653886e+02 *
        (2.54206000e+00 +
         T * (-2.75506200e-05 +
              T * (-3.10280300e-09 +
                   T * (4.55106700e-12 + T * -4.36805200e-16))));
    where(big_T, cp[index_v<OH>]) =
        4.88855960e+02 *
        (2.86472886e+00 +
         T * (1.05650448e-03 +
              T * (-2.59082758e-07 +
                   T * (3.05218674e-11 + T * -1.33195876e-15))));
    where(big_T, cp[index_v<H2>]) =
        4.12423719e+03 *
        (2.99142300e+00 +
         T * (7.00064400e-04 +
              T * (-5.63382900e-08 +
                   T * (-9.23157800e-12 + T * 1.58275200e-15))));
    where(big_T, cp[index_v<H2O>]) =
        4.61504339e+02 *
        (2.67214600e+00 +
         T * (3.05629300e-03 +
              T * (-8.73026000e-07 +
                   T * (1.20099600e-10 + T * -6.39161800e-15))));
    where(big_T, cp[index_v<HE>]) =
        2.07861554e+03 *
        (2.50000000e+00 +
         T * (0.00000000e+00 +
              T * (0.00000000e+00 +
                   T * (0.00000000e+00 + T * 0.00000000e+00))));

    where(big_T, cp[index_v<HO2>]) =
        2.51892334e+02 *
        (4.01721090e+00 +
         T * (2.23982013e-03 +
              T * (-6.33658150e-07 +
                   T * (1.14246370e-10 + T * -1.07908535e-14))));
    where(big_T, cp[index_v<H2O2>]) =
        2.44427980e+02 *
        (4.57316700e+00 +
         T * (4.33613600e-03 +
              T * (-1.47468900e-06 +
                   T * (2.34890400e-10 + T * -1.43165400e-14))));
  }
  auto not_big_T = !big_T;
  if (any_of(not_big_T)) {
    where(not_big_T, cp[index_v<N2>]) =
        2.96733125e+02 *
        (3.29867700e+00 +
         T * (1.40824000e-03 +
              T * (-3.96322200e-06 +
                   T * (5.64151500e-09 + T * -2.44485500e-12))));
    where(not_big_T, cp[index_v<Ar>]) =
        2.08132126e+02 *
        (2.50000000e+00 +
         T * (0.00000000e+00 +
              T * (0.00000000e+00 +
                   T * (0.00000000e+00 + T * 0.00000000e+00))));
    where(not_big_T, cp[index_v<H>]) =
        8.24847438e+03 *
        (2.50000000e+00 +
         T * (0.00000000e+00 +
              T * (0.00000000e+00 +
                   T * (0.00000000e+00 + T * 0.00000000e+00))));
    where(not_big_T, cp[index_v<O2>]) =
        2.59826943e+02 *
        (3.21293600e+00 +
         T * (1.12748600e-03 +
              T * (-5.75615000e-07 +
                   T * (1.31387700e-09 + T * -8.76855400e-13))));
    where(not_big_T, cp[index_v<O>]) =
        5.19653886e+02 *
        (2.94642900e+00 +
         T * (-1.63816600e-03 +
              T * (2.42103200e-06 +
                   T * (-1.60284300e-09 + T * 3.89069600e-13))));
    where(not_big_T, cp[index_v<OH>]) =
        4.88855960e+02 *
        (4.12530561e+00 +
         T * (-3.22544939e-03 +
              T * (6.52764691e-06 +
                   T * (-5.79853643e-09 + T * 2.06237379e-12))));
    where(not_big_T, cp[index_v<H2>]) =
        4.12423719e+03 *
        (3.29812400e+00 +
         T * (8.24944200e-04 +
              T * (-8.14301500e-07 +
                   T * (-9.47543400e-11 + T * 4.13487200e-13))));
    where(not_big_T, cp[index_v<H2O>]) =
        4.61504339e+02 *
        (3.38684200e+00 +
         T * (3.47498200e-03 +
              T * (-6.35469600e-06 +
                   T * (6.96858100e-09 + T * -2.50658800e-12))));
    where(not_big_T, cp[index_v<HE>]) =
        2.07861554e+03 *
        (2.50000000e+00 +
         T * (0.00000000e+00 +
              T * (0.00000000e+00 +
                   T * (0.00000000e+00 + T * 0.00000000e+00))));
    where(not_big_T, cp[index_v<HO2>]) =
        2.51892334e+02 *
        (4.30179801e+00 +
         T * (-4.74912051e-03 +
              T * (2.11582891e-05 +
                   T * (-2.42763894e-08 + T * 9.29225124e-12))));
    where(not_big_T, cp[index_v<H2O2>]) =
        2.44427980e+02 *
        (3.38875400e+00 +
         T * (6.56922600e-03 +
              T * (-1.48501300e-07 +
                   T * (-4.62580600e-09 + T * 2.47151500e-12))));
  }
}

template <typename Arithmetic>
std::enable_if_t<std::is_arithmetic<Arithmetic>::value>
get_specific_enthalpies_of_formation_impl(span<Arithmetic> h,
                                          Arithmetic T) noexcept {
  /*
          This function computes enthalpy 'h' and heat capacity 'cp' as
          function of temperature 'T' for all non steady state species
          in units [J/kg] and [J/kg K], respectively.
          The pArameter h and cp should provide workspace of length 11 */

  int i;
  if (T > 1000.0) {
    h[index_v<N2>] =
        2.96733125e+02 *
        (T * (2.92664000e+00 +
              T * (7.43988500e-04 +
                   T * (-1.89492033e-07 +
                        T * (2.52426000e-11 + T * -1.35067020e-15)))) -
         9.22797700e+02);
    h[index_v<Ar>] =
        2.08132126e+02 *
        (T * (2.50000000e+00 +
              T * (0.00000000e+00 +
                   T * (0.00000000e+00 +
                        T * (0.00000000e+00 + T * 0.00000000e+00)))) -
         7.45375000e+02);
    h[index_v<H>] =
        8.24847438e+03 *
        (T * (2.50000000e+00 +
              T * (0.00000000e+00 +
                   T * (0.00000000e+00 +
                        T * (0.00000000e+00 + T * 0.00000000e+00)))) +
         2.54716300e+04);
    h[index_v<O2>] =
        2.59826943e+02 *
        (T * (3.69757800e+00 +
              T * (3.06759850e-04 +
                   T * (-4.19614000e-08 +
                        T * (4.43820250e-12 + T * -2.27287000e-16)))) -
         1.23393000e+03);
    h[index_v<O>] =
        5.19653886e+02 *
        (T * (2.54206000e+00 +
              T * (-1.37753100e-05 +
                   T * (-1.03426767e-09 +
                        T * (1.13776675e-12 + T * -8.73610400e-17)))) +
         2.92308000e+04);
    h[index_v<OH>] =
        4.88855960e+02 *
        (T * (2.86472886e+00 +
              T * (5.28252240e-04 +
                   T * (-8.63609193e-08 +
                        T * (7.63046685e-12 + T * -2.66391752e-16)))) +
         3.68362875e+03);
    h[index_v<H2>] =
        4.12423719e+03 *
        (T * (2.99142300e+00 +
              T * (3.50032200e-04 +
                   T * (-1.87794300e-08 +
                        T * (-2.30789450e-12 + T * 3.16550400e-16)))) -
         8.35034000e+02);
    h[index_v<H2O>] =
        4.61504339e+02 *
        (T * (2.67214600e+00 +
              T * (1.52814650e-03 +
                   T * (-2.91008667e-07 +
                        T * (3.00249000e-11 + T * -1.27832360e-15)))) -
         2.98992100e+04);

    h[index_v<HE>] =
        2.07861554e+03 *
        (T * (2.50000000e+00 +
              T * (0.00000000e+00 +
                   T * (0.00000000e+00 +
                        T * (0.00000000e+00 + T * 0.00000000e+00)))) -
         7.45375000e+02);
    h[index_v<HO2>] =
        2.51892334e+02 *
        (T * (4.01721090e+00 +
              T * (1.11991006e-03 +
                   T * (-2.11219383e-07 +
                        T * (2.85615925e-11 + T * -2.15817070e-15)))) +
         1.11856713e+02);
    h[index_v<H2O2>] =
        2.44427980e+02 *
        (T * (4.57316700e+00 +
              T * (2.16806800e-03 +
                   T * (-4.91563000e-07 +
                        T * (5.87226000e-11 + T * -2.86330800e-15)))) -
         1.80069600e+04);
  } else if (T >= 299.999999) {
    h[index_v<N2>] =
        2.96733125e+02 *
        (T * (3.29867700e+00 +
              T * (7.04120000e-04 +
                   T * (-1.32107400e-06 +
                        T * (1.41037875e-09 + T * -4.88971000e-13)))) -
         1.02090000e+03);
    h[index_v<Ar>] =
        2.08132126e+02 *
        (T * (2.50000000e+00 +
              T * (0.00000000e+00 +
                   T * (0.00000000e+00 +
                        T * (0.00000000e+00 + T * 0.00000000e+00)))) -
         7.45375000e+02);
    h[index_v<H>] =
        8.24847438e+03 *
        (T * (2.50000000e+00 +
              T * (0.00000000e+00 +
                   T * (0.00000000e+00 +
                        T * (0.00000000e+00 + T * 0.00000000e+00)))) +
         2.54716300e+04);
    h[index_v<O2>] =
        2.59826943e+02 *
        (T * (3.21293600e+00 +
              T * (5.63743000e-04 +
                   T * (-1.91871667e-07 +
                        T * (3.28469250e-10 + T * -1.75371080e-13)))) -
         1.00524900e+03);

    h[index_v<O>] =
        5.19653886e+02 *
        (T * (2.94642900e+00 +
              T * (-8.19083000e-04 +
                   T * (8.07010667e-07 +
                        T * (-4.00710750e-10 + T * 7.78139200e-14)))) +
         2.91476400e+04);
    h[index_v<OH>] =
        4.88855960e+02 *
        (T * (4.12530561e+00 +
              T * (-1.61272470e-03 +
                   T * (2.17588230e-06 +
                        T * (-1.44963411e-09 + T * 4.12474758e-13)))) +
         3.34630913e+03);
    h[index_v<H2>] =
        4.12423719e+03 *
        (T * (3.29812400e+00 +
              T * (4.12472100e-04 +
                   T * (-2.71433833e-07 +
                        T * (-2.36885850e-11 + T * 8.26974400e-14)))) -
         1.01252100e+03);
    h[index_v<H2O>] =
        4.61504339e+02 *
        (T * (3.38684200e+00 +
              T * (1.73749100e-03 +
                   T * (-2.11823200e-06 +
                        T * (1.74214525e-09 + T * -5.01317600e-13)))) -
         3.02081100e+04);
    h[index_v<HE>] =
        2.07861554e+03 *
        (T * (2.50000000e+00 +
              T * (0.00000000e+00 +
                   T * (0.00000000e+00 +
                        T * (0.00000000e+00 + T * 0.00000000e+00)))) -
         7.45375000e+02);
    h[index_v<HO2>] =
        2.51892334e+02 *
        (T * (4.30179801e+00 +
              T * (-2.37456025e-03 +
                   T * (7.05276303e-06 +
                        T * (-6.06909735e-09 + T * 1.85845025e-12)))) +
         2.94808040e+02);
    h[index_v<H2O2>] =
        2.44427980e+02 *
        (T * (3.38875400e+00 +
              T * (3.28461300e-03 +
                   T * (-4.95004333e-08 +
                        T * (-1.15645150e-09 + T * 4.94303000e-13)))) -
         1.76631500e+04);
  } else {
    std::array<Arithmetic, 11> cp;
    get_specific_heat_capacities_at_constant_pressure_impl<Arithmetic>(cp,
                                                                       300.0);
    get_specific_enthalpies_of_formation_impl<Arithmetic>(h, 300.0);
    for (i = 0; i < h.size(); i++) {
      h[i] = (T - 300.) * cp[i] + h[i];
    }
  }
}

template <typename Simd>
std::enable_if_t<fub::is_simd<Simd>::value>
get_specific_enthalpies_of_formation_impl(span<Simd> h, Simd T) noexcept {

  auto big_T = (T > 1000.0);
  if (any_of(big_T)) {
    where(big_T, h[index_v<N2>]) =
        2.96733125e+02 *
        (T * (2.92664000e+00 +
              T * (7.43988500e-04 +
                   T * (-1.89492033e-07 +
                        T * (2.52426000e-11 + T * -1.35067020e-15)))) -
         9.22797700e+02);
    where(big_T, h[index_v<Ar>]) =
        2.08132126e+02 *
        (T * (2.50000000e+00 +
              T * (0.00000000e+00 +
                   T * (0.00000000e+00 +
                        T * (0.00000000e+00 + T * 0.00000000e+00)))) -
         7.45375000e+02);
    where(big_T, h[index_v<H>]) =
        8.24847438e+03 *
        (T * (2.50000000e+00 +
              T * (0.00000000e+00 +
                   T * (0.00000000e+00 +
                        T * (0.00000000e+00 + T * 0.00000000e+00)))) +
         2.54716300e+04);
    where(big_T, h[index_v<O2>]) =
        2.59826943e+02 *
        (T * (3.69757800e+00 +
              T * (3.06759850e-04 +
                   T * (-4.19614000e-08 +
                        T * (4.43820250e-12 + T * -2.27287000e-16)))) -
         1.23393000e+03);
    where(big_T, h[index_v<O>]) =
        5.19653886e+02 *
        (T * (2.54206000e+00 +
              T * (-1.37753100e-05 +
                   T * (-1.03426767e-09 +
                        T * (1.13776675e-12 + T * -8.73610400e-17)))) +
         2.92308000e+04);
    where(big_T, h[index_v<OH>]) =
        4.88855960e+02 *
        (T * (2.86472886e+00 +
              T * (5.28252240e-04 +
                   T * (-8.63609193e-08 +
                        T * (7.63046685e-12 + T * -2.66391752e-16)))) +
         3.68362875e+03);
    where(big_T, h[index_v<H2>]) =
        4.12423719e+03 *
        (T * (2.99142300e+00 +
              T * (3.50032200e-04 +
                   T * (-1.87794300e-08 +
                        T * (-2.30789450e-12 + T * 3.16550400e-16)))) -
         8.35034000e+02);
    where(big_T, h[index_v<H2O>]) =
        4.61504339e+02 *
        (T * (2.67214600e+00 +
              T * (1.52814650e-03 +
                   T * (-2.91008667e-07 +
                        T * (3.00249000e-11 + T * -1.27832360e-15)))) -
         2.98992100e+04);

    where(big_T, h[index_v<HE>]) =
        2.07861554e+03 *
        (T * (2.50000000e+00 +
              T * (0.00000000e+00 +
                   T * (0.00000000e+00 +
                        T * (0.00000000e+00 + T * 0.00000000e+00)))) -
         7.45375000e+02);
    where(big_T, h[index_v<HO2>]) =
        2.51892334e+02 *
        (T * (4.01721090e+00 +
              T * (1.11991006e-03 +
                   T * (-2.11219383e-07 +
                        T * (2.85615925e-11 + T * -2.15817070e-15)))) +
         1.11856713e+02);
    where(big_T, h[index_v<H2O2>]) =
        2.44427980e+02 *
        (T * (4.57316700e+00 +
              T * (2.16806800e-03 +
                   T * (-4.91563000e-07 +
                        T * (5.87226000e-11 + T * -2.86330800e-15)))) -
         1.80069600e+04);
  }

  auto mid_T = !big_T && (T >= 299.999999);
  if (any_of(mid_T)) {
    where(mid_T, h[index_v<N2>]) =
        2.96733125e+02 *
        (T * (3.29867700e+00 +
              T * (7.04120000e-04 +
                   T * (-1.32107400e-06 +
                        T * (1.41037875e-09 + T * -4.88971000e-13)))) -
         1.02090000e+03);
    where(mid_T, h[index_v<Ar>]) =
        2.08132126e+02 *
        (T * (2.50000000e+00 +
              T * (0.00000000e+00 +
                   T * (0.00000000e+00 +
                        T * (0.00000000e+00 + T * 0.00000000e+00)))) -
         7.45375000e+02);
    where(mid_T, h[index_v<H>]) =
        8.24847438e+03 *
        (T * (2.50000000e+00 +
              T * (0.00000000e+00 +
                   T * (0.00000000e+00 +
                        T * (0.00000000e+00 + T * 0.00000000e+00)))) +
         2.54716300e+04);
    where(mid_T, h[index_v<O2>]) =
        2.59826943e+02 *
        (T * (3.21293600e+00 +
              T * (5.63743000e-04 +
                   T * (-1.91871667e-07 +
                        T * (3.28469250e-10 + T * -1.75371080e-13)))) -
         1.00524900e+03);

    where(mid_T, h[index_v<O>]) =
        5.19653886e+02 *
        (T * (2.94642900e+00 +
              T * (-8.19083000e-04 +
                   T * (8.07010667e-07 +
                        T * (-4.00710750e-10 + T * 7.78139200e-14)))) +
         2.91476400e+04);
    where(mid_T, h[index_v<OH>]) =
        4.88855960e+02 *
        (T * (4.12530561e+00 +
              T * (-1.61272470e-03 +
                   T * (2.17588230e-06 +
                        T * (-1.44963411e-09 + T * 4.12474758e-13)))) +
         3.34630913e+03);
    where(mid_T, h[index_v<H2>]) =
        4.12423719e+03 *
        (T * (3.29812400e+00 +
              T * (4.12472100e-04 +
                   T * (-2.71433833e-07 +
                        T * (-2.36885850e-11 + T * 8.26974400e-14)))) -
         1.01252100e+03);
    where(mid_T, h[index_v<H2O>]) =
        4.61504339e+02 *
        (T * (3.38684200e+00 +
              T * (1.73749100e-03 +
                   T * (-2.11823200e-06 +
                        T * (1.74214525e-09 + T * -5.01317600e-13)))) -
         3.02081100e+04);
    where(mid_T, h[index_v<HE>]) =
        2.07861554e+03 *
        (T * (2.50000000e+00 +
              T * (0.00000000e+00 +
                   T * (0.00000000e+00 +
                        T * (0.00000000e+00 + T * 0.00000000e+00)))) -
         7.45375000e+02);
    where(mid_T, h[index_v<HO2>]) =
        2.51892334e+02 *
        (T * (4.30179801e+00 +
              T * (-2.37456025e-03 +
                   T * (7.05276303e-06 +
                        T * (-6.06909735e-09 + T * 1.85845025e-12)))) +
         2.94808040e+02);
    where(mid_T, h[index_v<H2O2>]) =
        2.44427980e+02 *
        (T * (3.38875400e+00 +
              T * (3.28461300e-03 +
                   T * (-4.95004333e-08 +
                        T * (-1.15645150e-09 + T * 4.94303000e-13)))) -
         1.76631500e+04);
  }

  auto small_T = !big_T && !mid_T;
  if (any_of(small_T)) {
    std::array<Simd, 11> cp;
    std::array<Simd, 11> h_;
    get_specific_heat_capacities_at_constant_pressure_impl<Simd>(cp, 300.0);
    get_specific_enthalpies_of_formation_impl<Simd>(h_, 300.0);
    for (int i = 0; i < 11; ++i) {
      where(small_T, h[i]) = (T - 300.) * cp[i] + h_[i];
    }
  }
}

} // namespace

/// @brief Fills the `cps` array with heat capacities which depend on
/// temperature [K] and mass fractions per each species.
void Burke2012::get_specific_heat_capacities_at_constant_pressure(
    species_span<double> cp, double temperature) noexcept {
  get_specific_heat_capacities_at_constant_pressure_impl<double>(cp,
                                                                 temperature);
}

/// @brief Fills the `cps` array with heat capacities which depend on
/// temperature [K] and mass fractions per each species.
void Burke2012::get_specific_heat_capacities_at_constant_pressure(
    simd_abi::scalar, species_span<scalar<double>> cp,
    scalar<double> temperature) noexcept {
  get_specific_heat_capacities_at_constant_pressure_impl<scalar<double>>(
      cp, temperature);
}

#if defined(FUB_SIMD_HAS_SSE)
void Burke2012::get_specific_heat_capacities_at_constant_pressure(
    simd_abi::sse, species_span<sse<double>> cp,
    sse<double> temperature) noexcept {
  get_specific_heat_capacities_at_constant_pressure_impl<sse<double>>(
      cp, temperature);
}
#endif

#if defined(FUB_SIMD_HAS_AVX)
void Burke2012::get_specific_heat_capacities_at_constant_pressure(
    simd_abi::avx, species_span<avx<double>> cp,
    avx<double> temperature) noexcept {
  get_specific_heat_capacities_at_constant_pressure_impl<avx<double>>(
      cp, temperature);
}
#endif

/// @brief Fills the `cps` array with heat capacities which depend on
/// temperature [K] and mass fractions per each species.
void Burke2012::get_specific_enthalpies_of_formation(
    species_span<double> enthalpies, double temperature) noexcept {
  get_specific_enthalpies_of_formation_impl<double>(enthalpies, temperature);
}

/// @brief Fills the `cps` array with heat capacities which depend on
/// temperature [K] and mass fractions per each species.
void Burke2012::get_specific_enthalpies_of_formation(
    simd_abi::scalar, species_span<scalar<double>> enthalpies,
    scalar<double> temperature) noexcept {
  get_specific_enthalpies_of_formation_impl<scalar<double>>(enthalpies,
                                                            temperature);
}

#if defined(FUB_SIMD_HAS_SSE)
/// @brief Fills the `cps` array with heat capacities which depend on
/// temperature [K] and mass fractions per each species.
void Burke2012::get_specific_enthalpies_of_formation(
    simd_abi::sse, species_span<sse<double>> enthalpies,
    sse<double> temperature) noexcept {
  get_specific_enthalpies_of_formation_impl<sse<double>>(enthalpies,
                                                         temperature);
}
#endif

#if defined(FUB_SIMD_HAS_AVX)
/// @brief Fills the `cps` array with heat capacities which depend on
/// temperature [K] and mass fractions per each species.
void Burke2012::get_specific_enthalpies_of_formation(
    simd_abi::avx, species_span<avx<double>> enthalpies,
    avx<double> temperature) noexcept {
  get_specific_enthalpies_of_formation_impl<avx<double>>(enthalpies,
                                                         temperature);
}
#endif

typedef enum SpeciesLabel {
  /* Computed species s.. */
  /* Steady-state species ss.. */
  sN2 = 0,
  sAR = 1,
  sH = 2,
  sO2 = 3,
  sO = 4,
  sOH = 5,
  sH2 = 6,
  sH2O = 7,
  sHE = 8,
  sHO2 = 9,
  sH2O2 = 10,
  sEnd
} SpeciesLabel;

typedef enum ReactionLabel {
  /* Reactions */
  r1f = 0,
  r1b = 1,
  r2f = 2,
  r2b = 3,
  r3f = 4,
  r3b = 5,
  r4f = 6,
  r4b = 7,
  r5f = 8,
  r5b = 9,
  r6f = 10,
  r6b = 11,
  r7f = 12,
  r7b = 13,
  r8f = 14,
  r8b = 15,
  r9f = 16,
  r9b = 17,
  r10f = 18,
  r10b = 19,
  r11f = 20,
  r11b = 21,
  r12f = 22,
  r12b = 23,
  r13f = 24,
  r13b = 25,
  r14f = 26,
  r14b = 27,
  r15f = 28,
  r15b = 29,
  r16f = 30,
  r16b = 31,
  r17f = 32,
  r17b = 33,
  r18f = 34,
  r18b = 35,
  r19f = 36,
  r19b = 37,
  r20f = 38,
  r20b = 39,
  r21f = 40,
  r21b = 41,
  r22f = 42,
  r22b = 43,
  r23f = 44,
  r23b = 45,
  r24f = 46,
  r24b = 47,
  r25f = 48,
  r25b = 49,
  r26f = 50,
  r26b = 51,
  r27f = 52,
  r27b = 53,
  /* PAHReactions */
  /* SootReactions */
  rEnd
} ReactionLabel;

typedef enum ThirdBodyLabel {
  mM1 = 0,
  mM2 = 1,
  mM3 = 2,
  mM4 = 3,
  mM5 = 4,
  mM6 = 5,
  mEnd
} ThirdBodyLabel;

namespace {
double GetLindRateCoeff(double temp, double pressure, double k0, double kInf,
                        double fc, double conc) {

  {
    const double R = 8314.34; /* [J / kmole K] */
    double Ntmp;
    double kl;
    double f;
    double cCoeff, dCoeff, log10kNull;

    int iTroe = 1;

    if (conc <= 0.0) {
      conc = pressure / (R * temp);
    }
    Ntmp = 0.75 - 1.27 * std::log10(fc);
    if (iTroe) {
      cCoeff = -0.4 - 0.67 * std::log10(fc);
      dCoeff = 0.14;
      k0 *= conc / std::max(kInf, 1.0e-60);
      log10kNull = std::log10(k0);
      f = (log10kNull + cCoeff) / (Ntmp - dCoeff * (log10kNull + cCoeff));
      f = std::pow(fc, 1.0 / (f * f + 1.0));
      kInf *= f * k0 / (1.0 + k0);
    } else {
      k0 = k0 * conc / kInf;
      kl = k0 / (1.0 + k0);
      f = std::log10(k0) / Ntmp;
      f = std::pow(fc, 1.0 / (f * f + 1.0));
      kInf = kInf * f * kl;
    }
    return kInf;
  }
}
} // namespace

void Burke2012::get_production_rates(span<double> cdot, span<const double> c,
                                     double temp, double pressure) noexcept {
  /*
     This function computes rates of production cdot in [kmole/(m^3s)].
     The parameters w ( reaction rate ), k ( rate coefficient )
     and M ( third body concentrations ) are just work space for this
     function. c contains the concentrations of non steady state species in
     [kmole/m^3] and is workspace for the steady state concentrations, which are
     computed in this function. temp is the temperature in [K] and pressure is
     the pressure in [Pa]. Called functions are 'GetLindRateCoeff',
     'ComputeSteadyStates', 'CatchZero' and the functions that evaluate the
     broadening factors of the Troe formulation of pressure dependant rate
     coefficients 'Fc*'
  */

  std::array<double, mEnd> M;
  std::array<double, rEnd> w;
  std::array<double, rEnd> k;

  // int nSpec = 11;
  // int nSpecIn = 11;
  double kTroe0, kTroeInf, fcTroe;
  double RGAS = 8314.34;
  double lgt = std::log(temp);
  double rt = RGAS * temp;

  M[mM1] = c[sN2] + c[sH] + c[sO2] + c[sO] + c[sOH] + 2.5 * c[sH2] +
           12 * c[sH2O] + c[sHO2] + c[sH2O2];

  M[mM2] = c[sN2] + c[sH] + c[sO2] + c[sO] + c[sOH] + 2.5 * c[sH2] +
           12 * c[sH2O] + c[sHO2] + c[sH2O2];

  M[mM3] = c[sN2] + 0.75 * c[sAR] + c[sH] + c[sO2] + c[sO] + c[sOH] +
           2.5 * c[sH2] + 12 * c[sH2O] + 0.75 * c[sHE] + c[sHO2] + c[sH2O2];

  M[mM4] = 2 * c[sN2] + c[sAR] + c[sH] + 1.5 * c[sO2] + c[sO] + c[sOH] +
           3 * c[sH2] + 1.1 * c[sHE] + c[sHO2] + c[sH2O2];

  M[mM5] = c[sN2] + 0.67 * c[sAR] + c[sH] + 0.78 * c[sO2] + c[sO] + c[sOH] +
           2 * c[sH2] + 14 * c[sH2O] + 0.8 * c[sHE] + c[sHO2] + c[sH2O2];

  M[mM6] = 1.5 * c[sN2] + c[sAR] + c[sH] + 1.2 * c[sO2] + c[sO] + c[sOH] +
           3.7 * c[sH2] + 7.5 * c[sH2O] + 0.65 * c[sHE] + c[sHO2] +
           7.7 * c[sH2O2];

  k[r1f] = 1.0400000000E+11 * std::exp(-63957000 / rt);
  k[r1b] = 2.0679845118E+08 * std::exp(0.434958 * lgt + 6605168.02 / rt);
  k[r2f] = 3.8180000000E+09 * std::exp(-33254000 / rt);
  k[r2b] = 2.2393817467E+09 * std::exp(-0.0350345 * lgt - 27487825.08 / rt);
  k[r3f] = 8.7920000000E+11 * std::exp(-80207000 / rt);
  k[r3b] = 5.1567952637E+11 * std::exp(-0.0350345 * lgt - 74440825.08 / rt);
  k[r4f] = 2.1600000000E+05 * std::exp(1.51 * lgt - 14351000 / rt);
  k[r4b] = 2.2182689717E+06 * std::exp(1.41638 * lgt - 77060225.89 / rt);
  k[r5f] = 3.3400000000E+01 * std::exp(2.42 * lgt + 8075000 / rt);
  k[r5b] = 5.8480989232E+02 * std::exp(2.36141 * lgt - 60400400.81 / rt);
  k[r6f] = 4.5770000000E+16 * std::exp(-1.4 * lgt - 436726000 / rt);
  k[r6b] = 1.9619373578E+14 * std::exp(-1.7488 * lgt - 3646272.681 / rt);
  k[r7f] = 5.8400000000E+15 * std::exp(-1.1 * lgt - 436726000 / rt);
  k[r7b] = 2.5033240484E+13 * std::exp(-1.4488 * lgt - 3646272.681 / rt);
  k[r8f] = 5.8400000000E+15 * std::exp(-1.1 * lgt - 436726000 / rt);
  k[r8b] = 2.5033240484E+13 * std::exp(-1.4488 * lgt - 3646272.681 / rt);
  k[r9f] = 6.1650000000E+09 * std::exp(-0.5 * lgt);
  k[r9b] = 4.2423561386E+14 * std::exp(-0.621192 * lgt - 497875720.4 / rt);
  k[r10f] = 1.8860000000E+07 * std::exp(7481000 / rt);
  k[r10b] = 1.2978237919E+12 * std::exp(-0.121192 * lgt - 490394720.4 / rt);
  k[r11f] = 1.8860000000E+07 * std::exp(7481000 / rt);
  k[r11b] = 1.2978237919E+12 * std::exp(-0.121192 * lgt - 490394720.4 / rt);
  k[r12f] = 4.7140000000E+12 * std::exp(-1 * lgt);
  k[r12b] = 6.4502650942E+14 * std::exp(-0.686234 * lgt - 427313552.4 / rt);
  k[r13f] = 6.0640000000E+24 * std::exp(-3.322 * lgt - 505385000 / rt);
  k[r13b] = 2.5310630492E+21 * std::exp(-3.57718 * lgt - 9596046.787 / rt);
  k[r14f] = 1.0060000000E+23 * std::exp(-2.44 * lgt - 502833000 / rt);
  k[r14b] = 4.1989601376E+19 * std::exp(-2.69518 * lgt - 7044046.787 / rt);
  kTroe0 = 6.3660000000E+14 * std::exp(-1.72 * lgt - 2196000 / rt);
  kTroeInf = 4.6510000000E+09 * std::exp(0.44 * lgt);
  fcTroe = 0.5 * std::exp(-temp / 1e-30) + 0.5 * std::exp(-temp / 1e+30);
  k[r15f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM5]);
  kTroe0 = 9.0292117899E+17 * std::exp(-1.73291 * lgt - 206862356.2 / rt);
  kTroeInf = 6.5967427010E+12 * std::exp(0.427093 * lgt - 204666356.2 / rt);
  k[r15b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM5]);
  k[r16f] = 2.7500000000E+03 * std::exp(2.09 * lgt + 6071000 / rt);
  k[r16b] = 4.5231945076E+02 * std::exp(2.45171 * lgt - 222342371.2 / rt);
  k[r17f] = 7.0790000000E+10 * std::exp(-1234000 / rt);
  k[r17b] = 1.3579714350E+07 * std::exp(0.761631 * lgt - 153319028.2 / rt);
  k[r18f] = 2.8500000000E+07 * std::exp(1 * lgt + 3029000 / rt);
  k[r18b] = 2.7494741434E+06 * std::exp(1.32667 * lgt - 219618196.2 / rt);
  k[r19f] = 2.8900000000E+10 * std::exp(2079000 / rt);
  k[r19b] = 4.8816975193E+10 * std::exp(0.268083 * lgt - 289043597.1 / rt);
  k[r20f] = 4.2000000000E+11 * std::exp(-50133000 / rt);
  k[r20b] = 6.3574616360E+13 * std::exp(-0.389273 * lgt - 212211053.7 / rt);
  k[r21f] = 1.3000000000E+08 * std::exp(6817000 / rt);
  k[r21b] = 1.9677857445E+10 * std::exp(-0.389273 * lgt - 155261053.7 / rt);
  kTroe0 = 2.4900000000E+21 * std::exp(-2.3 * lgt - 203970000 / rt);
  kTroeInf = 2.0000000000E+12 * std::exp(0.9 * lgt - 203966000 / rt);
  fcTroe = 0.57 * std::exp(-temp / 1e-30) + 0.43 * std::exp(-temp / 1e+30);
  k[r22f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM6]);
  kTroe0 = 2.2248483092E+12 * std::exp(-1.13619 * lgt + 10689381.66 / rt);
  kTroeInf = 1.7870267543E+03 * std::exp(2.06381 * lgt + 10693381.66 / rt);
  k[r22b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM6]);
  k[r23f] = 2.4100000000E+10 * std::exp(-16610000 / rt);
  k[r23b] = 5.1591045676E+04 * std::exp(1.41899 * lgt - 297739571.6 / rt);
  k[r24f] = 4.8200000000E+10 * std::exp(-33263000 / rt);
  k[r24b] = 5.2375134409E+07 * std::exp(0.75098 * lgt - 99598317.44 / rt);
  k[r25f] = 9.5500000000E+03 * std::exp(2 * lgt - 16610000 / rt);
  k[r25b] = 6.0865850326E+00 * std::exp(2.71595 * lgt - 77179142.52 / rt);
  k[r26f] = 1.7400000000E+09 * std::exp(-1331000 / rt);
  k[r26b] = 1.9417254097E+07 * std::exp(0.657356 * lgt - 130375543.3 / rt);
  k[r27f] = 7.5900000000E+10 * std::exp(-30418000 / rt);
  k[r27b] = 8.4699401492E+08 * std::exp(0.657356 * lgt - 159462543.3 / rt);

  w[r1f] = k[r1f] * c[sH] * c[sO2];
  w[r1b] = k[r1b] * c[sOH] * c[sO];
  w[r2f] = k[r2f] * c[sO] * c[sH2];
  w[r2b] = k[r2b] * c[sOH] * c[sH];
  w[r3f] = k[r3f] * c[sO] * c[sH2];
  w[r3b] = k[r3b] * c[sOH] * c[sH];
  w[r4f] = k[r4f] * c[sH2] * c[sOH];
  w[r4b] = k[r4b] * c[sH] * c[sH2O];
  w[r5f] = k[r5f] * c[sOH] * c[sOH];
  w[r5b] = k[r5b] * c[sH2O] * c[sO];
  w[r6f] = k[r6f] * c[sH2] * M[mM1];
  w[r6b] = k[r6b] * c[sH] * c[sH] * M[mM1];
  w[r7f] = k[r7f] * c[sH2] * c[sAR];
  w[r7b] = k[r7b] * c[sAR] * c[sH] * c[sH];
  w[r8f] = k[r8f] * c[sH2] * c[sHE];
  w[r8b] = k[r8b] * c[sHE] * c[sH] * c[sH];
  w[r9f] = k[r9f] * c[sO] * c[sO] * M[mM2];
  w[r9b] = k[r9b] * c[sO2] * M[mM2];
  w[r10f] = k[r10f] * c[sO] * c[sO] * c[sAR];
  w[r10b] = k[r10b] * c[sAR] * c[sO2];
  w[r11f] = k[r11f] * c[sO] * c[sO] * c[sHE];
  w[r11b] = k[r11b] * c[sHE] * c[sO2];
  w[r12f] = k[r12f] * c[sO] * c[sH] * M[mM3];
  w[r12b] = k[r12b] * c[sOH] * M[mM3];
  w[r13f] = k[r13f] * c[sH2O] * M[mM4];
  w[r13b] = k[r13b] * c[sOH] * c[sH] * M[mM4];
  w[r14f] = k[r14f] * c[sH2O] * c[sH2O];
  w[r14b] = k[r14b] * c[sH2O] * c[sOH] * c[sH];
  w[r15f] = k[r15f] * c[sH] * c[sO2];
  w[r15b] = k[r15b] * c[sHO2];
  w[r16f] = k[r16f] * c[sHO2] * c[sH];
  w[r16b] = k[r16b] * c[sO2] * c[sH2];
  w[r17f] = k[r17f] * c[sHO2] * c[sH];
  w[r17b] = k[r17b] * c[sOH] * c[sOH];
  w[r18f] = k[r18f] * c[sHO2] * c[sO];
  w[r18b] = k[r18b] * c[sOH] * c[sO2];
  w[r19f] = k[r19f] * c[sHO2] * c[sOH];
  w[r19b] = k[r19b] * c[sO2] * c[sH2O];
  w[r20f] = k[r20f] * c[sHO2] * c[sHO2];
  w[r20b] = k[r20b] * c[sO2] * c[sH2O2];
  w[r21f] = k[r21f] * c[sHO2] * c[sHO2];
  w[r21b] = k[r21b] * c[sO2] * c[sH2O2];
  w[r22f] = k[r22f] * c[sH2O2];
  w[r22b] = k[r22b] * c[sOH] * c[sOH];
  w[r23f] = k[r23f] * c[sH2O2] * c[sH];
  w[r23b] = k[r23b] * c[sOH] * c[sH2O];
  w[r24f] = k[r24f] * c[sH2O2] * c[sH];
  w[r24b] = k[r24b] * c[sH2] * c[sHO2];
  w[r25f] = k[r25f] * c[sH2O2] * c[sO];
  w[r25b] = k[r25b] * c[sHO2] * c[sOH];
  w[r26f] = k[r26f] * c[sH2O2] * c[sOH];
  w[r26b] = k[r26b] * c[sH2O] * c[sHO2];
  w[r27f] = k[r27f] * c[sH2O2] * c[sOH];
  w[r27b] = k[r27b] * c[sH2O] * c[sHO2];

  cdot[sN2] = 0.0;

  cdot[sAR] = -w[r7f] + w[r7f] - w[r7b] + w[r7b] - w[r10f] + w[r10f] - w[r10b] +
              w[r10b];

  cdot[sH] = -w[r1f] + w[r1b] + w[r2f] - w[r2b] + w[r3f] - w[r3b] + w[r4f] -
             w[r4b] + 2 * w[r6f] - 2 * w[r6b] + 2 * w[r7f] - 2 * w[r7b] +
             2 * w[r8f] - 2 * w[r8b] - w[r12f] + w[r12b] + w[r13f] - w[r13b] +
             w[r14f] - w[r14b] - w[r15f] + w[r15b] - w[r16f] + w[r16b] -
             w[r17f] + w[r17b] - w[r23f] + w[r23b] - w[r24f] + w[r24b];

  cdot[sO2] = -w[r1f] + w[r1b] + w[r9f] - w[r9b] + w[r10f] - w[r10b] + w[r11f] -
              w[r11b] - w[r15f] + w[r15b] + w[r16f] - w[r16b] + w[r18f] -
              w[r18b] + w[r19f] - w[r19b] + w[r20f] - w[r20b] + w[r21f] -
              w[r21b];

  cdot[sO] = w[r1f] - w[r1b] - w[r2f] + w[r2b] - w[r3f] + w[r3b] + w[r5f] -
             w[r5b] - 2 * w[r9f] + 2 * w[r9b] - 2 * w[r10f] + 2 * w[r10b] -
             2 * w[r11f] + 2 * w[r11b] - w[r12f] + w[r12b] - w[r18f] + w[r18b] -
             w[r25f] + w[r25b];

  cdot[sOH] = w[r1f] - w[r1b] + w[r2f] - w[r2b] + w[r3f] - w[r3b] - w[r4f] +
              w[r4b] - 2 * w[r5f] + 2 * w[r5b] + w[r12f] - w[r12b] + w[r13f] -
              w[r13b] + w[r14f] - w[r14b] + 2 * w[r17f] - 2 * w[r17b] +
              w[r18f] - w[r18b] - w[r19f] + w[r19b] + 2 * w[r22f] -
              2 * w[r22b] + w[r23f] - w[r23b] + w[r25f] - w[r25b] - w[r26f] +
              w[r26b] - w[r27f] + w[r27b];

  cdot[sH2] = -w[r2f] + w[r2b] - w[r3f] + w[r3b] - w[r4f] + w[r4b] - w[r6f] +
              w[r6b] - w[r7f] + w[r7b] - w[r8f] + w[r8b] + w[r16f] - w[r16b] +
              w[r24f] - w[r24b];

  cdot[sH2O] = w[r4f] - w[r4b] + w[r5f] - w[r5b] - w[r13f] + w[r13b] -
               2 * w[r14f] + w[r14f] - w[r14b] + 2 * w[r14b] + w[r19f] -
               w[r19b] + w[r23f] - w[r23b] + w[r26f] - w[r26b] + w[r27f] -
               w[r27b];

  cdot[sHE] = -w[r8f] + w[r8f] - w[r8b] + w[r8b] - w[r11f] + w[r11f] - w[r11b] +
              w[r11b];

  cdot[sHO2] = w[r15f] - w[r15b] - w[r16f] + w[r16b] - w[r17f] + w[r17b] -
               w[r18f] + w[r18b] - w[r19f] + w[r19b] - 2 * w[r20f] +
               2 * w[r20b] - 2 * w[r21f] + 2 * w[r21b] + w[r24f] - w[r24b] +
               w[r25f] - w[r25b] + w[r26f] - w[r26b] + w[r27f] - w[r27b];

  cdot[sH2O2] = w[r20f] - w[r20b] + w[r21f] - w[r21b] - w[r22f] + w[r22b] -
                w[r23f] + w[r23b] - w[r24f] + w[r24b] - w[r25f] + w[r25b] -
                w[r26f] + w[r26b] - w[r27f] + w[r27b];
}

} // namespace burke2012
} // namespace mechanism
} // namespace euler
} // namespace fub
