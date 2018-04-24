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

#ifndef FUB_EULER_MECHANISM_GRI_30_HPP
#define FUB_EULER_MECHANISM_GRI_30_HPP

#include "fub/euler/ideal_gas.hpp"

#include <tuple>

namespace fub {
namespace euler {
namespace mechanism {
namespace gri_30 {

////////////////////////////////////////////////////////////////////////////////
// List of Species
// {{{

struct AR {
  static constexpr const char* name() noexcept { return "AR"; }
};

struct O {
  static constexpr const char* name() noexcept { return "O"; }
};

struct O2 {
  static constexpr const char* name() noexcept { return "O2"; }
};

struct H {
  static constexpr const char* name() noexcept { return "H"; }
};

struct OH {
  static constexpr const char* name() noexcept { return "OH"; }
};

struct H2 {
  static constexpr const char* name() noexcept { return "H2"; }
};

struct HO2 {
  static constexpr const char* name() noexcept { return "HO2"; }
};

struct H2O2 {
  static constexpr const char* name() noexcept { return "H2O2"; }
};

struct CH {
  static constexpr const char* name() noexcept { return "CH"; }
};

struct CO {
  static constexpr const char* name() noexcept { return "CO"; }
};

struct _3XCH2 {
  static constexpr const char* name() noexcept { return "3XCH2"; }
};

struct HCO {
  static constexpr const char* name() noexcept { return "HCO"; }
};

struct _1XCH2 {
  static constexpr const char* name() noexcept { return "1XCH2"; }
};

struct CH3 {
  static constexpr const char* name() noexcept { return "CH3"; }
};

struct CH2O {
  static constexpr const char* name() noexcept { return "CH2O"; }
};

struct CH4 {
  static constexpr const char* name() noexcept { return "CH4"; }
};

struct CO2 {
  static constexpr const char* name() noexcept { return "CO2"; }
};

struct CH2OH {
  static constexpr const char* name() noexcept { return "CH2OH"; }
};

struct CH3O {
  static constexpr const char* name() noexcept { return "CH3O"; }
};

struct CH3OH {
  static constexpr const char* name() noexcept { return "CH3OH"; }
};

struct C2H {
  static constexpr const char* name() noexcept { return "C2H"; }
};

struct C2H2 {
  static constexpr const char* name() noexcept { return "C2H2"; }
};

struct HCCO {
  static constexpr const char* name() noexcept { return "HCCO"; }
};

struct C2H3 {
  static constexpr const char* name() noexcept { return "C2H3"; }
};

struct CH2CO {
  static constexpr const char* name() noexcept { return "CH2CO"; }
};

struct C2H4 {
  static constexpr const char* name() noexcept { return "C2H4"; }
};

struct C2H5 {
  static constexpr const char* name() noexcept { return "C2H5"; }
};

struct C2H6 {
  static constexpr const char* name() noexcept { return "C2H6"; }
};

struct H2O {
  static constexpr const char* name() noexcept { return "H2O"; }
};

struct N2 {
  static constexpr const char* name() noexcept { return "N2"; }
};

struct C {
  static constexpr const char* name() noexcept { return "C"; }
};

struct HCCOH {
  static constexpr const char* name() noexcept { return "HCCOH"; }
};

struct N {
  static constexpr const char* name() noexcept { return "N"; }
};

struct NO {
  static constexpr const char* name() noexcept { return "NO"; }
};

struct N2O {
  static constexpr const char* name() noexcept { return "N2O"; }
};

struct NO2 {
  static constexpr const char* name() noexcept { return "NO2"; }
};

struct NH {
  static constexpr const char* name() noexcept { return "NH"; }
};

struct HNO {
  static constexpr const char* name() noexcept { return "HNO"; }
};

struct NH2 {
  static constexpr const char* name() noexcept { return "NH2"; }
};

struct NNH {
  static constexpr const char* name() noexcept { return "NNH"; }
};

struct CN {
  static constexpr const char* name() noexcept { return "CN"; }
};

struct NCO {
  static constexpr const char* name() noexcept { return "NCO"; }
};

struct HCN {
  static constexpr const char* name() noexcept { return "HCN"; }
};

struct HOCN {
  static constexpr const char* name() noexcept { return "HOCN"; }
};

struct HNCO {
  static constexpr const char* name() noexcept { return "HNCO"; }
};

struct H2CN {
  static constexpr const char* name() noexcept { return "H2CN"; }
};

struct HCNN {
  static constexpr const char* name() noexcept { return "HCNN"; }
};

struct HCNO {
  static constexpr const char* name() noexcept { return "HCNO"; }
};

struct NH3 {
  static constexpr const char* name() noexcept { return "NH3"; }
};

struct CH2CHO {
  static constexpr const char* name() noexcept { return "CH2CHO"; }
};

struct CH3CHO {
  static constexpr const char* name() noexcept { return "CH3CHO"; }
};

struct C3H8 {
  static constexpr const char* name() noexcept { return "C3H8"; }
};

struct C3H7 {
  static constexpr const char* name() noexcept { return "C3H7"; }
};
// }}}

////////////////////////////////////////////////////////////////////////////////
//                                                                [class.gri_30]
struct gri_30 {
  using species_tuple =
      std::tuple<AR, O, O2, H, OH, H2, HO2, H2O2, CH, CO, _3XCH2, HCO, _1XCH2,
                 CH3, CH2O, CH4, CO2, CH2OH, CH3O, CH3OH, C2H, C2H2, HCCO, C2H3,
                 CH2CO, C2H4, C2H5, C2H6, H2O, N2, C, HCCOH, N, NO, N2O, NO2,
                 NH, HNO, NH2, NNH, CN, NCO, HCN, HOCN, HNCO, H2CN, HCNN, HCNO,
                 NH3, CH2CHO, CH3CHO, C3H8, C3H7>;

  template <typename T>
  using species_span = span<T, std::tuple_size<species_tuple>::value>;

  /// @brief Returns the mass per mole in [g/mol].
  ///
  /// These are constant values measured from experiments.
  static species_span<const double> get_molar_masses() noexcept;

  /// @brief Returns dY/dt for a specified mass distribution, temperature and
  /// pressure.
  static void get_production_rates(species_span<double> dydt,
                                   species_span<const double> y,
                                   double temperature,
                                   double pressure) noexcept;

  /// @brief Fills the `cps` array with heat capacities which depend on
  /// temperature [K] and mass fractions per each species.
  static void get_specific_heat_capacities_at_constant_pressure(
      species_span<double> cp, double temperature) noexcept;

  /// @brief Fills the `cps` array with heat capacities which depend on
  /// temperature [K] and mass fractions per each species.
  static void get_specific_heat_capacities_at_constant_pressure(
      simd_abi::scalar, species_span<scalar<double>> cp,
      scalar<double> temperature) noexcept;

#if defined(FUB_SIMD_HAS_SSE)
  /// @brief Fills the `cps` array with heat capacities which depend on
  /// temperature [K] and mass fractions per each species.
  static void get_specific_heat_capacities_at_constant_pressure(
      simd_abi::sse, species_span<sse<double>> cp,
      sse<double> temperature) noexcept;
#endif

#if defined(FUB_SIMD_HAS_AVX)
  /// @brief Fills the `cps` array with heat capacities which depend on
  /// temperature [K] and mass fractions per each species.
  static void get_specific_heat_capacities_at_constant_pressure(
      simd_abi::avx, species_span<avx<double>> cp,
      avx<double> temperature) noexcept;
#endif

  /// @brief Fills the `cps` array with heat capacities which depend on
  /// temperature [K] and mass fractions per each species.
  static void
  get_specific_enthalpies_of_formation(species_span<double> enthalpies,
                                       double temperature) noexcept;

  /// @brief Fills the `cps` array with heat capacities which depend on
  /// temperature [K] and mass fractions per each species.
  static void
  get_specific_enthalpies_of_formation(simd_abi::scalar,
                                       species_span<scalar<double>> enthalpies,
                                       scalar<double> temperature) noexcept;

#if defined(FUB_SIMD_HAS_SSE)
  /// @brief Fills the `cps` array with heat capacities which depend on
  /// temperature [K] and mass fractions per each species.
  static void
  get_specific_enthalpies_of_formation(simd_abi::sse,
                                       species_span<sse<double>> enthalpies,
                                       sse<double> temperature) noexcept;
#endif

#if defined(FUB_SIMD_HAS_AVX)
  /// @brief Fills the `cps` array with heat capacities which depend on
  /// temperature [K] and mass fractions per each species.
  static void
  get_specific_enthalpies_of_formation(simd_abi::avx,
                                       species_span<avx<double>> enthalpies,
                                       avx<double> temperature) noexcept;
#endif
};

template <typename T>
static constexpr int Index_v =
    variable_find_index<T, gri_30::species_tuple>::value;

template <typename T> constexpr int as_index(T) noexcept { return Index_v<T>; }

} // namespace gri_30
} // namespace mechanism
} // namespace euler
} // namespace fub

#endif // !BURKE2012_HPP
