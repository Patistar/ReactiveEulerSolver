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

#ifndef FUB_EULER_MECHANISM_BURKE2012_HPP
#define FUB_EULER_MECHANISM_BURKE2012_HPP

#include "fub/euler/ideal_gas.hpp"

#include <tuple>

namespace fub {
namespace euler {
namespace mechanism {
namespace burke2012 {

////////////////////////////////////////////////////////////////////////////////
// List of Species
// {{{
struct N2 {
  using value_type = double;
  static constexpr const char* name() noexcept { return "N2"; }
};
struct Ar {
  using value_type = double;
  static constexpr const char* name() noexcept { return "Ar"; }
};
struct H {
  using value_type = double;
  static constexpr const char* name() noexcept { return "H"; }
};
struct O2 {
  using value_type = double;
  static constexpr const char* name() noexcept { return "O2"; }
};
struct O {
  using value_type = double;
  static constexpr const char* name() noexcept { return "O"; }
};
struct OH {
  using value_type = double;
  static constexpr const char* name() noexcept { return "OH"; }
};
struct H2 {
  using value_type = double;
  static constexpr const char* name() noexcept { return "H2"; }
};
struct H2O {
  using value_type = double;
  static constexpr const char* name() noexcept { return "H2O"; }
};
struct HE {
  using value_type = double;
  static constexpr const char* name() noexcept { return "HE"; }
};
struct HO2 {
  using value_type = double;
  static constexpr const char* name() noexcept { return "HO2"; }
};
struct H2O2 {
  using value_type = double;
  static constexpr const char* name() noexcept { return "H2O2"; }
};

namespace variables {
static constexpr N2 n2{};
static constexpr Ar ar{};
static constexpr H h{};
static constexpr O2 o2{};
static constexpr O o{};
static constexpr OH oh{};
static constexpr H2 h2{};
static constexpr H2O h2o{};
static constexpr HE he{};
static constexpr HO2 ho2{};
static constexpr H2O2 h2o2{};
} // namespace variables

// }}}

////////////////////////////////////////////////////////////////////////////////
//                                                             [class.Burke2012]
struct Burke2012 {
  using species_tuple =
      std::tuple<N2, Ar, H, O2, O, OH, H2, H2O, HE, HO2, H2O2>;

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
    variable_find_index<T, Burke2012::species_tuple>::value;

template <typename T> constexpr int as_index(T) noexcept { return Index_v<T>; }

} // namespace burke2012
} // namespace mechanism
} // namespace euler
} // namespace fub

#endif // !BURKE2012_HPP
