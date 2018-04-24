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

#ifndef FUB_EULER_MECHANISM_SINGLE_STAGE_HPP
#define FUB_EULER_MECHANISM_SINGLE_STAGE_HPP

/// @file This file containts a simple reaction mechanism for 3 species and 2
/// reactions.
///
/// The details are explained in the thesis of P. Berndt.

#include "fub/euler/ideal_gas.hpp"
#include "fub/tuple.hpp"
#include "fub/utility.hpp"

namespace fub {
namespace euler {
namespace mechanism {
namespace single_stage {
struct Fuel {
  using value_type = double;
  static constexpr const char* name() noexcept { return "Fuel"; }
};
struct Intermediate {
  using value_type = double;
  static constexpr const char* name() noexcept { return "Intermediate"; }
};
struct Product {
  using value_type = double;
  static constexpr const char* name() noexcept { return "Product"; }
};

namespace variables {
static constexpr Fuel fuel{};
static constexpr Intermediate intermediate{};
static constexpr Product product{};
} // namespace variables

struct single_stage {
  ///@brief Defines all species.
  using species_tuple = std::tuple<Fuel, Intermediate, Product>;
  static constexpr int species_size = std::tuple_size<species_tuple>::value;
  template <typename T> using species_span = span<T, species_size>;
};

template <typename T>
static constexpr int Index_v =
    variable_find_index<T, single_stage::species_tuple>::value;

template <typename T> constexpr int as_index(T) noexcept { return Index_v<T>; }

} // namespace single_stage
} // namespace mechanism

template <int Rank>
class ideal_gas<mechanism::single_stage::single_stage, Rank, double> {
public:
  //////////////////////////////////////////////////////////////////////////////
  // Typedefs
  // {{{
  /// We use this placeholder to determine the precision in which we do
  /// computations.
  /// @todo Make this configurable.
  using arithmetic_type = double;
  using mechanism_type = mechanism::single_stage::single_stage;

  /// @brief This is a variables type which contains all species.
  ///
  /// It is of the form
  ///
  ///    `Species = std::tuple<Species0, Species1, ..., SpeciesN>`
  ///
  /// and thus also contains the Species type which is only implicitly given
  /// in the State variables.
  using species_tuple = std::tuple<mechanism::single_stage::Fuel,
                                   mechanism::single_stage::Intermediate,
                                   mechanism::single_stage::Product>;
  static constexpr std::size_t species_size =
      std::tuple_size<species_tuple>::value;

  template <typename T> using species_span = span<T, species_size>;

  struct momentum_states {
    template <typename> struct impl;
    template <int... Is> struct impl<std::integer_sequence<int, Is...>> {
      using type = quantities<Momentum<Is>...>;
    };
    using type = typename impl<std::make_integer_sequence<int, Rank>>::type;
  };
  using momentum_states_t = typename momentum_states::type;

  /// @brief This defines all variables on a data patch for this problem.
  /// Note that it does not contain the first species concentration.
  /// It is given implicitly by `ρY_0 = ρ(1 - Σ Y_i)`
  using complete_state =
      concat_t<quantities<Density>, momentum_states_t,
               quantities<Energy, Pressure, Temperature, SpeedOfSound>,
               tail_t<as_quantities_t<species_tuple>>>;

  /// @brief This defines all variables which shall be conserved by the finite
  /// volume methods.
  using conservative_state =
      concat_t<quantities<Density>, momentum_states_t, quantities<Energy>,
               tail_t<as_quantities_t<species_tuple>>>;

  /// @brief We assume a constant heat capacity ratio of 1.28
  static constexpr double gamma = 1.28;
  /// @brief Heat capacity at constant volume, this follows from gamma and R
  static constexpr double cv = 1. / (gamma - 1.);
  /// @brief Heat capacity at constant pressure, this follows from gamma and R
  static constexpr double cp = 1 + cv;
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // get molar masses
  // {{{
  /// @brief Returns an array filled with Ones {1, 1, 1}.
  static species_span<const double> get_molar_masses() noexcept {
    static constexpr std::array<double, species_size> m{{1.0, 1.0, 1.0}};
    return m;
  }
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // get specific enthalpies_of_formation
  // {{{
  /// @brief Fills the `enthalpies` array with zeros.
  /// @{
  static void
  get_specific_enthalpies_of_formation(species_span<double> enthalpies,
                                       double) noexcept {
    std::fill(enthalpies.begin(), enthalpies.end(), 0.0);
  }

  template <typename Abi>
  static void get_specific_enthalpies_of_formation(
      const Abi&, species_span<simd<double, Abi>> enthalpies,
      simd<double, Abi>) noexcept {
    std::fill(enthalpies.begin(), enthalpies.end(), simd<double, Abi>(0.0));
  }
  /// @}
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // get specific heat capacity at constant pressure
  // {{{

  /// @brief Fills the `cps` array with the constant `cp`.
  /// @{
  static void
  get_specific_heat_capacities_at_constant_pressure(species_span<double> cps,
                                                    double) noexcept {
    std::fill(cps.begin(), cps.end(), cp);
  }

  template <typename Abi>
  static std::enable_if_t<is_simd_abi_v<Abi>>
  get_specific_heat_capacities_at_constant_pressure(
      Abi, species_span<simd<double, Abi>> cps, simd<double, Abi>) noexcept {
    std::fill(cps.begin(), cps.end(), cp);
  }
  /// @}
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // get reference temperature
  static constexpr double get_reference_temperature() noexcept { return 0.0; }

  //////////////////////////////////////////////////////////////////////////////
  // get universal gas constant
  static constexpr double get_universal_gas_constant() noexcept { return 1.0; }

  //////////////////////////////////////////////////////////////////////////////
  // get reference specific enthalpies of formation
  // {{{
  static constexpr double get_mean_reference_specific_enthalpies_of_formation(
      species_span<const double>) noexcept {
    return 0;
  }
  template <typename Abi>
  std::enable_if_t<is_simd_abi_v<Abi>, simd<double, Abi>>
  get_mean_reference_specific_enthalpy_of_formation(
      const Abi&, nodeduce_t<species_span<const simd<double, Abi>>>) const
      noexcept {
    return 0;
  }
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // get specific gas constants
  // {{{
  static species_span<const double> get_specific_gas_constants() noexcept {
    static constexpr std::array<double, 3> gas_constants{{1., 1., 1.}};
    return gas_constants;
  }
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // get specific internal energy
  // {{{

  /// @brief Computes the internal energy dependent on temperature.
  ///
  /// We are using the formula
  ///
  ///     U = \int c_v(T) dT = cv T
  /// @{
  double get_specific_internal_energy(double temperature) const noexcept {
    return cv * temperature;
  }

  template <typename Abi>
  std::enable_if_t<is_simd_abi_v<Abi>, simd<double, Abi>>
  get_specific_internal_energy(const Abi&, simd<double, Abi> temperature) const
      noexcept {
    return cv * temperature;
  }
  /// @}
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // set URY
  // {{{
  /// @brief Returns a valid and complete set of variables computed from
  /// conservative variables alone.
  ///
  /// This function will be used to retain a complete state after integrating
  /// the conservative variables in time
  /// @{
  complete_state set_URY(double internal_energy, double rho,
                         species_span<const double> mass_fractions) const {
    using namespace variables;
    complete_state q;
    q[density] = rho;
    q[temperature] = internal_energy / cv;
    q[energy] = rho * internal_energy;
    q[pressure] = rho * q[temperature];
    q[speed_of_sound] = std::sqrt(gamma * q[temperature]);
    fub::for_each_tuple_element(
        [&](auto species) {
          q[species] = rho * mass_fractions[as_index(species)];
        },
        tail_t<species_tuple>{});
    return q;
  }

  template <typename Abi>
  std::enable_if_t<is_simd_abi_v<Abi>, add_simd_t<complete_state, Abi>> set_URY(
      const Abi&, simd<double, Abi> internal_energy, simd<double, Abi> rho,
      nodeduce_t<species_span<const simd<double, Abi>>> mass_fractions) const {
    using namespace variables;
    add_simd_t<complete_state, Abi> q;
    q[density] = rho;
    q[temperature] = internal_energy / cv;
    q[energy] = rho * internal_energy;
    q[pressure] = rho * q[temperature];
    q[speed_of_sound] = fub::sqrt(gamma * q[temperature]);
    fub::for_each_tuple_element(
        [&](auto species) {
          q[species] = rho * mass_fractions[as_index(species)];
        },
        tail_t<species_tuple>{});
    return q;
  }
  // }}}
  /// @}

  //////////////////////////////////////////////////////////////////////////////
  // set TPY / TPX
  // {{{
  complete_state set_TPY(double T, double p, species_span<const double> Y) const
      noexcept {
    assert(T >= 0);
    assert(p >= 0);
    complete_state q;
    using namespace variables;
    q[temperature] = T;
    q[pressure] = p;
    q[density] = p / T;
    q[energy] = p * cv;
    q[speed_of_sound] = std::sqrt(gamma * T);
    fub::for_each_tuple_element(
        [&](auto species) { q[species] = q[density] * Y[as_index(species)]; },
        tail_t<species_tuple>());
    return q;
  }

  complete_state set_TPY(species_span<double> hs, species_span<double> cps,
                         double T, double p, species_span<const double> Y) const
      noexcept {

    using namespace variables;
    get_specific_enthalpies_of_formation(hs, T);
    get_specific_heat_capacities_at_constant_pressure(cps, T);
    return set_TPY(T, p, Y);
  }

  template <typename Abi>
  std::enable_if_t<is_simd_abi_v<Abi>, add_simd_t<complete_state, Abi>>
  set_TPY(const Abi&, simd<double, Abi> T, simd<double, Abi> p,
          nodeduce_t<species_span<const simd<double, Abi>>> Y) const noexcept {
    add_simd_t<complete_state, Abi> q;
    using namespace variables;
    q[temperature] = T;
    q[pressure] = p;
    q[density] = p / T;
    q[energy] = p * cv;
    q[speed_of_sound] = fub::sqrt(gamma * T);
    fub::for_each_tuple_element(
        [&](auto species) { q[species] = q[density] * Y[as_index(species)]; },
        tail_t<species_tuple>());
    return q;
  }

  template <typename Abi>
  std::enable_if_t<is_simd_abi_v<Abi>, add_simd_t<complete_state, Abi>>
  set_TPY(const Abi& abi, species_span<simd<double, Abi>> hs,
          species_span<simd<double, Abi>> cps, simd<double, Abi> T,
          simd<double, Abi> p, species_span<const simd<double, Abi>> Y) const
      noexcept {
    using namespace variables;
    get_specific_enthalpies_of_formation(abi, hs, T);
    get_specific_heat_capacities_at_constant_pressure(abi, cps, T);
    return set_TPY(abi, T, p, Y);
  }
  std::array<double, species_size>
  get_mass_from_mole_fractions(species_span<const double> mole_fractions) const
      noexcept {
    std::array<double, species_size> fractions;
    species_span<const double> molar_masses = get_molar_masses();
    std::transform(mole_fractions.begin(), mole_fractions.end(),
                   molar_masses.begin(), fractions.begin(),
                   std::multiplies<>{});
    double sum = std::accumulate(fractions.begin(), fractions.end(), double(0));
    std::transform(fractions.begin(), fractions.end(), fractions.begin(),
                   [sum](double Y) { return Y / sum; });
    return fractions;
  }

  complete_state set_TPX(species_span<double> hs, species_span<double> cps,
                         double T, double p,
                         species_span<const double> mole_fractions) const
      noexcept {
    std::array<double, species_size> mass_fractions =
        get_mass_from_mole_fractions(mole_fractions);
    return set_TPY(hs, cps, T, p, mass_fractions);
  }

  complete_state set_TPX(double T, double p,
                         species_span<const double> mole_fractions) const
      noexcept {
    std::array<double, species_size> mass_fractions =
        get_mass_from_mole_fractions(mole_fractions);
    return set_TPY(T, p, mass_fractions);
  }

  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // set momentum / velocity
  // {{{
  template <typename State>
  std::array<double, Rank> get_momentum(const State& q) const noexcept {
    return fub::apply(
        [&q](auto... dims) {
          return std::array<double, Rank>{
              {q[variables::momentum<decltype(dims)::value>]...}};
        },
        as_tuple_t<std::make_integer_sequence<int, Rank>>());
  }

  template <typename Abi, typename State>
  std::array<simd<double, Abi>, Rank> get_momentum(const Abi&,
                                                   const State& q) const
      noexcept {
    return fub::apply(
        [&q](auto... dims) {
          return std::array<simd<double, Abi>, Rank>{
              {q[variables::momentum<decltype(dims)::value>]...}};
        },
        as_tuple_t<std::make_integer_sequence<int, Rank>>());
  }

  template <typename State>
  State set_momentum(const State& q, const std::array<double, Rank>& rhou) const
      noexcept {
    using namespace variables;
    State result{q};
    std::array<double, Rank> old = get_momentum(q);
    auto e_kinetic_new = fub::transform_reduce(rhou, rhou, double(0));
    auto e_kinetic_old = fub::transform_reduce(old, old, double(0));
    result[energy] += 0.5 * (e_kinetic_new - e_kinetic_old) / q[density];
    for_each_tuple_element(
        [&](auto dir) { result[momentum<decltype(dir)::value>] = rhou[dir()]; },
        as_tuple_t<std::make_integer_sequence<int, Rank>>());
    return result;
  }

  template <typename Abi, typename State>
  std::enable_if_t<is_simd_abi_v<Abi>, State> set_momentum(
      const Abi& abi, const State& q,
      const nodeduce_t<std::array<simd<double, Abi>, Rank>>& rhou) const
      noexcept {
    using namespace variables;
    State result{q};
    std::array<simd<double, Abi>, Rank> old = get_momentum(abi, q);
    auto e_kinetic_new =
        fub::transform_reduce(rhou, rhou, simd<double, Abi>(0));
    auto e_kinetic_old = fub::transform_reduce(old, old, simd<double, Abi>(0));
    result[energy] += 0.5 * (e_kinetic_new - e_kinetic_old) / q[density];
    for_each_tuple_element(
        [&](auto dir) { result[momentum<decltype(dir)::value>] = rhou[dir()]; },
        as_tuple_t<std::make_integer_sequence<int, Rank>>());
    return result;
  }

  template <typename State>
  State set_velocity(const State& q, std::array<double, Rank> u) const
      noexcept {
    using variables::density;
    std::transform(u.begin(), u.end(), u.begin(),
                   [rho = q[density]](double u_i) { return rho * u_i; });
    return set_momentum(q, u);
  }

  template <typename Abi, typename State>
  std::enable_if_t<is_simd_abi_v<Abi>, State>
  set_velocity(const Abi& abi, const State& q,
               nodeduce_t<std::array<simd<double, Abi>, Rank>> u) const
      noexcept {
    using variables::density;
    std::transform(
        u.begin(), u.end(), u.begin(),
        [rho = q[density]](simd<double, Abi> u_i) { return rho * u_i; });
    return set_momentum(abi, q, u);
  }
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // get flux
  // {{{

  /// @brief Returns and evaluates the euler flux function for conservative
  /// variables for a specified complete state `q`.
  ///
  /// This function simply advects all quantities but energy and momentum,
  /// which additionally get their respective pressure contributions.
  template <int Dim>
  add_flux_t<conservative_state> get_flux(int_constant<Dim>,
                                          const complete_state& q) const
      noexcept {
    using namespace variables;
    add_flux_t<conservative_state> f;
    const auto u = q[momentum<Dim>] / q[density];
    fub::for_each_tuple_element([&](auto fq) { f[fq] = u * q[unflux(fq)]; },
                                get_variables(f));
    f[flux(momentum<Dim>)] += q[pressure];
    f[flux(energy)] += u * q[pressure];
    return f;
  }

  add_flux_t<conservative_state> get_flux(const complete_state& q) const
      noexcept {
    return get_flux(int_c<0>, q);
  }

  template <int Dim, typename Abi>
  add_flux_t<add_simd_t<conservative_state, Abi>>
  get_flux(int_constant<Dim>, const Abi&,
           const add_simd_t<complete_state, Abi>& q) const noexcept {
    using namespace variables;
    add_flux_t<add_simd_t<conservative_state, Abi>> f;
    auto u = q[momentum<Dim>] / q[density];
    fub::for_each_tuple_element([&](auto fq) { f[fq] = u * q[unflux(fq)]; },
                                get_variables(f));
    f[flux(momentum<Dim>)] += q[pressure];
    f[flux(energy)] += u * q[pressure];
    return f;
  }

  template <typename Abi>
  add_flux_t<add_simd_t<conservative_state, Abi>>
  get_flux(const Abi& abi, const add_simd_t<complete_state, Abi>& q) const
      noexcept {
    return get_flux(int_c<0>, abi, q);
  }
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // get mass fractions from a state
  // {{{
  template <typename State>
  static std::array<double, species_size>
  get_mass_fractions(const State& q) noexcept {
    std::array<double, species_size> mass_fractions;
    fub::for_each_tuple_element(
        [&](auto species) {
          using variables::density;
          mass_fractions[as_index(species)] = q[species] / q[density];
        },
        tail_t<species_tuple>());
    mass_fractions[as_index(head_t<species_tuple>())] = std::max(
        double(0), 1 - std::accumulate(mass_fractions.begin() + 1,
                                       mass_fractions.end(), double(0)));
    return mass_fractions;
  }

  template <typename Abi, typename State>
  static std::array<simd<double, Abi>, species_size>
  get_mass_fractions(const Abi&, const State& q) noexcept {
    std::array<simd<double, Abi>, species_size> mass_fractions;
    fub::for_each_tuple_element(
        [&](auto species) {
          using variables::density;
          mass_fractions[as_index(species)] = q[species] / q[density];
        },
        tail_t<species_tuple>());
    mass_fractions[as_index(head_t<species_tuple>())] = fub::max(
        simd<double, Abi>(0),
        1 - std::accumulate(mass_fractions.begin() + 1, mass_fractions.end(),
                            simd<double, Abi>(0)));
    return mass_fractions;
  }
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // get mean specific enthalpy of formation from a state
  // {{{
  template <typename State>
  double get_specific_enthalpies_of_formation(const State&) const noexcept {
    return 0;
  }
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // get mean specific heat capacity at constant pressure from a state
  // {{{
  double get_mean_specific_heat_capacity_at_constant_pressure(
      const complete_state&) const noexcept {
    return cp;
  }

  double get_mean_specific_heat_capacity_at_constant_pressure(
      const conservative_state&) const noexcept {
    return cp;
  }

  template <typename Abi>
  simd<double, Abi> get_mean_specific_heat_capacity_at_constant_pressure(
      Abi, const add_simd_t<complete_state, Abi>&) const noexcept {
    return cp;
  }
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // get a complete state from a conservative one
  // {{{

  /// @brief Retains a complete state from a conservative one.
  ///
  /// This function maybe used after a state update of conservative variables.
  complete_state get_complete_state(const conservative_state& q) const {
    using namespace variables;
    conservative_state Q = set_momentum(q, {0});
    std::array<double, species_size> mass_fractions = get_mass_fractions(Q);
    double U = Q[energy] / Q[density];
    auto thermo = set_URY(U, Q[density], mass_fractions);
    return set_momentum(thermo, get_momentum(q));
  }

  template <typename Abi>
  add_simd_t<complete_state, Abi>
  get_complete_state(const Abi& abi,
                     const add_simd_t<conservative_state, Abi>& q) const {
    using namespace variables;
    add_simd_t<conservative_state, Abi> Q = set_momentum(abi, q, {});
    auto mass_fractions = get_mass_fractions(abi, Q);
    simd<double, Abi> U = Q[energy] / Q[density];
    add_simd_t<complete_state, Abi> thermo =
        set_URY(abi, U, Q[density], mass_fractions);
    return set_momentum(abi, thermo, get_momentum(abi, q));
  }
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // getter for variables
  // {{{
  template <typename Var>
  double get(Var, const complete_state& q) const noexcept {
    return q.template get<Var>();
  }

  template <int Dir>
  double get(Velocity<Dir>, const complete_state& q) const noexcept {
    using namespace variables;
    return q[momentum<Dir>] / q[density];
  }

  double get(head_t<species_tuple>, const complete_state& q) const noexcept {
    using namespace variables;
    double implicit_fraction =
        q[density] -
        foldl(tail_t<species_tuple>(), double(0),
              [&](double sum, auto species) { return sum + q[species]; });
    return std::max(double(0), implicit_fraction);
  }

  template <typename Abi, typename Var>
  simd<double, Abi> get(Abi, Var,
                        const add_simd_t<complete_state, Abi>& q) const
      noexcept {
    return q.template get<Var>();
  }

  template <int Dir, typename Abi>
  simd<double, Abi> get(Abi, Velocity<Dir>,
                        const add_simd_t<complete_state, Abi>& q) const
      noexcept {
    using namespace variables;
    return q[momentum<Dir>] / q[density];
  }

  template <typename Abi>
  simd<double, Abi> get(Abi, head_t<species_tuple>,
                        const add_simd_t<complete_state, Abi>& q) const
      noexcept {
    using namespace variables;
    simd<double, Abi> implicit_fraction =
        q[density] - foldl(tail_t<species_tuple>(), simd<double, Abi>(0),
                           [&](simd<double, Abi> sum, auto species) {
                             return sum + q[species];
                           });
    return fub::max(simd<double, Abi>(0), implicit_fraction);
  }
  // }}}
};

/// @brief We assume a constant heat capacity ratio of 1.28
template <int Rank>
constexpr double
    ideal_gas<mechanism::single_stage::single_stage, Rank, double>::gamma;
/// @brief Heat capacity at constant volume, this follows from gamma and R
template <int Rank>
constexpr double
    ideal_gas<mechanism::single_stage::single_stage, Rank, double>::cv;
/// @brief Heat capacity at constant pressure, this follows from gamma and R
template <int Rank>
constexpr double
    ideal_gas<mechanism::single_stage::single_stage, Rank, double>::cp;

} // namespace euler
} // namespace fub

#endif // !SINGLESTATE_HPP
