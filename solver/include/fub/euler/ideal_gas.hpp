// Copyright(c) 2018 Maikel Nadolski
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

#ifndef FUB_EULER_IDEAL_GAS_HPP
#define FUB_EULER_IDEAL_GAS_HPP

#include "fub/euler/variables.hpp"
#include "fub/simd.hpp"
#include "fub/span.hpp"
#include "fub/tuple.hpp"
#include "fub/variables.hpp"

#include <numeric>
#include <tuple>

namespace fub {
namespace euler {

/// @brief This problem describes the multiple species ideal gas euler
/// equations.
///
/// @details The system is given by s + 3 equations for s + 1 species
///
///    ρ_t    + (ρu)_x         = 0
///    (ρu)_t + (ρu^2 + p)_x   = 0
///    (ρE)_t + (u (ρE + p))_x = 0
///    (ρY)_t + (ρuY)_x        = 0           (Y = (Y_1, ... Y_s))
///
/// and we consider the ideal gas equation of state
///
///    p = ρRT
///
/// The first species mass fraction is implicitly given by `Y_0 = 1 - Σ Y_i`.
/// We parametrize this equation of state by a reaction mechanism which can be
/// used via the \a MechanismTraits class to do all computations.
///
/// @note This is not necessarily a equation of state for perfect gases and
/// heat capacities may depend on temperature.
template <typename Mechanism, int Rank = 1, typename A = double>
class ideal_gas {
private:
#if __has_cpp_attribute(no_unique_address)
  [[no_unique_address]] Mechanism m_mechanism;
#else
  Mechanism m_mechanism;
#endif

public:
  //////////////////////////////////////////////////////////////////////////////
  // Typedefs
  // {{{
  /// We use this placeholder to determine the precision in which we do
  /// computations.
  /// @todo Make this configurable.
  using arithmetic_type = A;
  static_assert(std::is_arithmetic<A>::value, "A is not an arithmetic type.");

  /// @brief This is a variables type which contains all species.
  ///
  /// It is of the form
  ///
  ///    `Species = std::tuple<Species0, Species1, ..., SpeciesN>`
  ///
  /// and thus also contains the Species type which is only implicitly given
  /// in the State variables.
  using species_tuple = typename Mechanism::species_tuple;
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
               quantities<Energy, Pressure, Temperature, SpeedOfSound,
                          HeatCapacityAtConstantP, Enthalpy>,
               tail_t<as_quantities_t<species_tuple>>>;

  /// @brief This defines all variables which shall be conserved by the finite
  /// volume methods.
  using conservative_state =
      concat_t<quantities<Density>, momentum_states_t, quantities<Energy>,
               tail_t<as_quantities_t<species_tuple>>>;
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // get molar masses
  // {{{
  species_span<const A> get_molar_masses() const noexcept {
    return m_mechanism.get_molar_masses();
  }
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // get production rates
  // {{{
  void get_production_rates(species_span<A> dYdt, species_span<const A> Y, A T,
                            A P) const noexcept {
    m_mechanism.get_production_rates(dYdt, Y, T, P);
  }
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // get specific enthalpies of formation
  // {{{
  void get_specific_enthalpies_of_formation(species_span<A> enthalpies,
                                            A temperature) const noexcept {
    m_mechanism.get_specific_enthalpies_of_formation(enthalpies, temperature);
  }

  template <typename Abi>
  std::enable_if_t<is_simd_abi_v<Abi>> get_specific_enthalpies_of_formation(
      const Abi& abi, nodeduce_t<species_span<simd<A, Abi>>> enthalpies,
      nodeduce_t<simd<A, Abi>> temperature) const noexcept {
    m_mechanism.get_specific_enthalpies_of_formation(abi, enthalpies,
                                                     temperature);
  }
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // get specific heat capacity at constant pressure
  // {{{
  void get_specific_heat_capacities_at_constant_pressure(
      species_span<A> heat_capacities, A temperature) const noexcept {
    assert(temperature > 0);
    m_mechanism.get_specific_heat_capacities_at_constant_pressure(
        heat_capacities, temperature);
  }

  template <typename Abi>
  std::enable_if_t<is_simd_abi_v<Abi>>
  get_specific_heat_capacities_at_constant_pressure(
      const Abi& abi, nodeduce_t<species_span<simd<A, Abi>>> heat_capacities,
      nodeduce_t<simd<A, Abi>> temperature) const noexcept {
    m_mechanism.get_specific_heat_capacities_at_constant_pressure(
        abi, heat_capacities, temperature);
  }
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // get specific gas constants
  // {{{
private:
  /// @brief Returns an array which contains all specific gas constants.
  auto store_specific_gas_constants() const noexcept {
    // [J / (kmol K)]
    const A universal_gas_constant = get_universal_gas_constant();
    // [kg / kmol]
    species_span<const A> molar_masses = get_molar_masses();

    std::array<A, species_size> constants;
    std::transform(molar_masses.begin(), molar_masses.end(), constants.begin(),
                   [=](A molar_mass) {
                     // [J / (kg K)]
                     return universal_gas_constant / molar_mass;
                   });
    return constants;
  }

public:
  /// @brief Returns a view to a static array containing specific gas constants
  /// for each species.
  species_span<const A> get_specific_gas_constants() const noexcept {
    static const auto specific_gas_constants = store_specific_gas_constants();
    return specific_gas_constants;
  }
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // get reference temperature
  // {{{
  constexpr A get_reference_temperature() const noexcept { return A(300); }
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // get reference specific enthalpies of formation
  // {{{
private:
  auto store_reference_specific_enthalpies_of_formation() const noexcept {
    const A T_ref = get_reference_temperature();

    std::array<A, species_size> constants;
    get_specific_enthalpies_of_formation(constants, T_ref);
    return constants;
  }

public:
  /// @brief Returns the specific enthalpy at reference temeprature for the
  /// specified mass fractions.
  A get_mean_reference_specific_enthalpy_of_formation(
      species_span<const A> mass_fractions) const noexcept {
    static const auto h0 = store_reference_specific_enthalpies_of_formation();
    return fub::transform_reduce(h0, mass_fractions, A(0));
  }

  /// @brief Returns the specific enthalpy at reference temperature for the
  /// specified mass fractions.
  template <typename Abi>
  std::enable_if_t<is_simd_abi_v<Abi>, simd<A, Abi>>
  get_mean_reference_specific_enthalpy_of_formation(
      const Abi&,
      nodeduce_t<species_span<const simd<A, Abi>>> mass_fractions) const
      noexcept {
    static const auto h0 = store_reference_specific_enthalpies_of_formation();
    return fub::transform_reduce(h0, mass_fractions, simd<A, Abi>(0));
  }
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // get universal gas constant
  // {{{
  A get_universal_gas_constant() const noexcept { return A(8314.45); }
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // get specific internal energy
  // {{{

  /// @brief Computes the internal energy dependent on temperature, enthalpy
  /// and mass fractions.
  ///
  /// We are using the formula
  ///
  ///     U = \int c_v(T) dT = \int (c_p(T) - R) dT = H - RT.
  ///
  /// Where H is the mean enthalpy, R the specific gas constant and `T`
  /// temperature.
  A get_specific_internal_energy(A temperature, A mean_specific_gas_constant,
                                 species_span<const A> enthalpies,
                                 species_span<const A> mass_fractions) const
      noexcept {
    // [J / kg]
    A mean_H = fub::transform_reduce(enthalpies, mass_fractions, A(0));
    return mean_H - mean_specific_gas_constant * temperature;
  }

  /// @brief Computes the internal energy dependent on temperature, enthalpy
  /// and mass fractions.
  ///
  /// We are using the formula
  ///
  ///     U = \int c_v(T) dT = \int (c_p(T) - R) dT = H - RT.
  ///
  /// Where H is the mean enthalpy, R the specific gas constant and `T`
  /// temperature.
  template <typename Abi>
  std::enable_if_t<is_simd_abi_v<Abi>, simd<A, Abi>>
  get_specific_internal_energy(
      const Abi&, nodeduce_t<simd<A, Abi>> temperature,
      nodeduce_t<simd<A, Abi>> mean_specific_gas_constant,
      nodeduce_t<species_span<const simd<A, Abi>>> enthalpies,
      nodeduce_t<species_span<const simd<A, Abi>>> mass_fractions) const
      noexcept {
    // [J / kg]
    simd<A, Abi> mean_H =
        fub::transform_reduce(enthalpies, mass_fractions, simd<A, Abi>(0));
    return mean_H - mean_specific_gas_constant * temperature;
  }
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // set URY
  // {{{
private:
  /// @brief Returns the temperature which matches a specified intenral
  /// energy.
  ///
  /// This function uses a Newton-like iteration to find a temperature such
  /// that
  ///
  ///     get_internal_energy(T, Y_0) == U_0
  ///
  /// for specified Y_0 and U_0.
  A find_temperature(species_span<A> h, species_span<A> cps, A U_0,
                     species_span<const A> Y_0, A T_0 = 300,
                     std::ptrdiff_t tolerance = 100000,
                     std::ptrdiff_t max_iterations = 1000) const {
    const species_span<const A> Rs = get_specific_gas_constants();
    const A R = fub::transform_reduce(Rs, Y_0, A(0));
    auto get_U = [&](A T) {
      get_specific_enthalpies_of_formation(h, T);
      get_specific_heat_capacities_at_constant_pressure(cps, T);
      return get_specific_internal_energy(T, R, h, Y_0);
    };
    A T = T_0;
    A U = get_U(T);
    std::ptrdiff_t iterations = 0;
    while (!almost_equal(U_0, U, tolerance) && iterations < max_iterations) {
      const A cp = fub::transform_reduce(cps, Y_0, A(0));
      const A cv = cp - R;
      const A dU = U_0 - U;
      const A dT = fub::clamp(dU / cv, A(-100), A(100));
      T = T + dT;
      U = get_U(T);
      ++iterations;
    }
    if (almost_equal(U_0, U, tolerance)) {
      return T;
    }
    throw std::runtime_error{
        "ideal_gas::find_temperature: Newton iteration did not converge."};
  }

  template <typename Abi>
  std::enable_if_t<is_simd_abi_v<Abi>, simd<A, Abi>>
  find_temperature(const Abi& abi, nodeduce_t<species_span<simd<A, Abi>>> h,
                   nodeduce_t<species_span<simd<A, Abi>>> cps, simd<A, Abi> U_0,
                   species_span<const simd<A, Abi>> Y_0,
                   simd<A, Abi> T_0 = 300., std::ptrdiff_t tolerance = 100000,
                   std::ptrdiff_t max_iterations = 1000) const {
    // Prepare U evaluation
    const species_span<const A> Rs = get_specific_gas_constants();
    const simd<A, Abi> R = fub::transform_reduce(Rs, Y_0, simd<A, Abi>(0));
    auto get_U = [&](simd<A, Abi> T) {
      get_specific_enthalpies_of_formation(abi, h, T);
      get_specific_heat_capacities_at_constant_pressure(abi, cps, T);
      return get_specific_internal_energy(abi, T, R, h, Y_0);
    };

    // Newton Iteration
    simd<A, Abi> T = T_0;
    simd<A, Abi> U = get_U(T);
    std::ptrdiff_t iterations = 0;
    bool not_done_yet = any_of(!almost_equal(U_0, U, tolerance));
    while (not_done_yet && iterations < max_iterations) {
      const simd<A, Abi> cp = fub::transform_reduce(cps, Y_0, simd<A, Abi>(0));
      const simd<A, Abi> cv = cp - R;
      const simd<A, Abi> dU = U_0 - U;
      const simd<A, Abi> dT = fub::clamp(dU / cv, A(-100), A(100));
      T = T + dT;
      U = get_U(T);
      ++iterations;
      not_done_yet = any_of(!almost_equal(U_0, U, tolerance));
    }
    // If done, return.
    if (!not_done_yet) {
      return T;
    }
    throw std::runtime_error{"ideal_gas::find_temperature (simd): Newton "
                             "iteration did not converge."};
  }

public:
  /// @brief Returns a valid and complete set of variables computed from
  /// conservative variables alone.
  ///
  /// This function will be used to retain a complete state after integrating
  /// the conservative variables in time
  complete_state set_URY(A internal_energy, A rho,
                         species_span<const A> mass_fractions) const {
    using namespace variables;
    complete_state q;
    q[density] = rho;
    array<A, species_size> hs;
    array<A, species_size> cps;
    q[temperature] = find_temperature(hs, cps, internal_energy, mass_fractions);
    q[cp] = fub::transform_reduce(cps, mass_fractions, A(0));
    q[enthalpy] = fub::transform_reduce(hs, mass_fractions, A(0));
    species_span<const A> specific_gas_constants = get_specific_gas_constants();
    const A R =
        fub::transform_reduce(specific_gas_constants, mass_fractions, A(0));
    q[pressure] = rho * R * q[temperature];
    const A gamma = q[cp] / (q[cp] - R);
    q[speed_of_sound] = std::sqrt(gamma * R * q[temperature]);
    fub::for_each_tuple_element(
        [&](auto species) {
          q[species] = rho * mass_fractions[as_index(species)];
        },
        tail_t<species_tuple>{});
    const A h_ref =
        get_mean_reference_specific_enthalpy_of_formation(mass_fractions);
    const A T_ref = get_reference_temperature();
    q[energy] = rho * (internal_energy + R * T_ref - h_ref);
    return q;
  }

  /// @brief Returns a valid and complete set of variables computed from
  /// conservative variables alone.
  ///
  /// This function will be used to retain a complete state after integrating
  /// the conservative variables in time
  template <typename Abi>
  std::enable_if_t<is_simd_abi_v<Abi>, add_simd_t<complete_state, Abi>>
  set_URY(const Abi& abi, nodeduce_t<simd<A, Abi>> internal_energy,
          nodeduce_t<simd<A, Abi>> rho,
          nodeduce_t<species_span<const simd<A, Abi>>> mass_fractions) const {
    using namespace variables;
    add_simd_t<complete_state, Abi> q{};
    q[density] = rho;
    array<simd<A, Abi>, species_size> hs;
    array<simd<A, Abi>, species_size> cps;
    simd<A, Abi> temp =
        find_temperature(abi, hs, cps, internal_energy, mass_fractions);
    q[temperature] = temp;
    q[enthalpy] = fub::transform_reduce(hs, mass_fractions, simd<A, Abi>(0));
    q[cp] = fub::transform_reduce(cps, mass_fractions, simd<A, Abi>(0));
    const species_span<const A> Rs = get_specific_gas_constants();
    const simd<A, Abi> gas_constant =
        fub::transform_reduce(Rs, mass_fractions, simd<A, Abi>(0));
    q[pressure] = rho * gas_constant * q[temperature];
    const simd<A, Abi> gamma = q[cp] / (q[cp] - gas_constant);
    q[speed_of_sound] = fub::sqrt(gamma * gas_constant * q[temperature]);
    fub::for_each_tuple_element(
        [&](auto species) {
          q[species] = rho * mass_fractions[as_index(species)];
        },
        tail_t<species_tuple>{});
    const simd<A, Abi> h_ref =
        get_mean_reference_specific_enthalpy_of_formation(abi, mass_fractions);
    const simd<A, Abi> T_ref = get_reference_temperature();
    q[energy] = rho * (internal_energy + gas_constant * T_ref - h_ref);
    return q;
  }

  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // set TPY / TPX
  // {{{
  std::array<A, species_size>
  get_mass_from_mole_fractions(species_span<const A> mole_fractions) const
      noexcept {
    std::array<A, species_size> fractions;
    species_span<const A> molar_masses = get_molar_masses();
    std::transform(mole_fractions.begin(), mole_fractions.end(),
                   molar_masses.begin(), fractions.begin(),
                   std::multiplies<>{});
    A sum = std::accumulate(fractions.begin(), fractions.end(), A(0));
    std::transform(fractions.begin(), fractions.end(), fractions.begin(),
                   [sum](A Y) { return Y / sum; });
    return fractions;
  }

  complete_state set_TPY(species_span<A> hs, species_span<A> cps, A T, A p,
                         species_span<const A> mass_fractions) const noexcept {
    using namespace variables;
    complete_state q;
    q[temperature] = T;
    q[pressure] = p;
    species_span<const A> Rs = get_specific_gas_constants();
    const A R = fub::transform_reduce(Rs, mass_fractions, A(0));
    q[density] = p / (R * T);
    {
      // Set total energy
      get_specific_enthalpies_of_formation(hs, T);
      q[enthalpy] = fub::transform_reduce(hs, mass_fractions, A(0));
      const A U = get_specific_internal_energy(T, R, hs, mass_fractions);
      const A T_ref = get_reference_temperature();
      const A h_ref =
          get_mean_reference_specific_enthalpy_of_formation(mass_fractions);
      const A E_ref = R * T_ref - h_ref;
      q[energy] = q[density] * (U + E_ref);
    }
    get_specific_heat_capacities_at_constant_pressure(cps, T);
    q[cp] = fub::transform_reduce(cps, mass_fractions, A(0));
    const A gamma = q[cp] / (q[cp] - R);
    assert(T >= 0);
    assert(gamma >= 0);
    assert(R >= 0);
    q[speed_of_sound] = std::sqrt(gamma * R * T);
    fub::for_each_tuple_element(
        [&q, &mass_fractions](auto species) {
          q[species] = q[density] * mass_fractions[as_index(species)];
        },
        tail_t<species_tuple>{});
    return q;
  }

  complete_state set_TPY(A T, A p, species_span<const A> mass_fractions) const
      noexcept {
    array<A, species_size> hs;
    array<A, species_size> cps;
    return set_TPY(hs, cps, T, p, mass_fractions);
  }

  template <typename Abi>
  add_simd_t<complete_state, Abi> set_TPY(
      const Abi& abi, const nodeduce_t<simd<A, Abi>>& T,
      const nodeduce_t<simd<A, Abi>>& p,
      const nodeduce_t<species_span<const simd<A, Abi>>>& mass_fractions) const
      noexcept {
    using namespace variables;
    add_simd_t<complete_state, Abi> q{};
    q[temperature] = T;
    q[pressure] = p;
    species_span<const A> Rs = get_specific_gas_constants();
    const simd<A, Abi> R =
        fub::transform_reduce(Rs, mass_fractions, simd<A, Abi>(0));
    q[density] = p / (R * T);
    {
      // Set total energy
      array<simd<A, Abi>, species_size> h;
      get_specific_enthalpies_of_formation(abi, h, T);
      q[enthalpy] = fub::transform_reduce(h, mass_fractions, simd<A, Abi>(0.0));
      const simd<A, Abi> U =
          get_specific_internal_energy(abi, T, R, h, mass_fractions);
      const simd<A, Abi> h_ref =
          get_mean_reference_specific_enthalpy_of_formation(abi,
                                                            mass_fractions);
      const A T_ref = get_reference_temperature();
      const simd<A, Abi> E_ref = R * T_ref - h_ref;
      q[energy] = q[density] * (U + E_ref);
    }
    array<simd<A, Abi>, species_size> cps;
    get_specific_heat_capacities_at_constant_pressure(abi, cps, T);
    q[cp] = fub::transform_reduce(cps, mass_fractions, simd<A, Abi>(0.0));
    const simd<A, Abi> gamma = q[cp] / (q[cp] - R);
    q[speed_of_sound] = fub::sqrt(gamma * R * T);
    fub::for_each_tuple_element(
        [&q, &mass_fractions](auto species) {
          q[species] = q[density] * mass_fractions[as_index(species)];
        },
        tail_t<species_tuple>{});
    return q;
  }

  complete_state set_TPX(species_span<A> hs, species_span<A> cps, A T, A p,
                         species_span<const A> mole_fractions) const noexcept {
    std::array<A, species_size> mass_fractions =
        get_mass_from_mole_fractions(mole_fractions);
    return set_TPY(hs, cps, T, p, mass_fractions);
  }

  complete_state set_TPX(A T, A p, species_span<const A> mole_fractions) const
      noexcept {

    std::array<A, species_size> mass_fractions =
        get_mass_from_mole_fractions(mole_fractions);
    return set_TPY(T, p, mass_fractions);
  }
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // set momentum / velocity
  // {{{
  template <typename State>
  std::array<A, Rank> get_momentum(const State& q) const noexcept {
    return fub::apply(
        [&q](auto... dims) {
          return std::array<A, Rank>{
              {q[variables::momentum<decltype(dims)::value>]...}};
        },
        as_tuple_t<std::make_integer_sequence<int, Rank>>());
  }

  template <typename Abi, typename State>
  std::array<simd<A, Abi>, Rank> get_momentum(const Abi&, const State& q) const
      noexcept {
    return fub::apply(
        [&q](auto... dims) {
          return std::array<simd<A, Abi>, Rank>{
              {q[variables::momentum<decltype(dims)::value>]...}};
        },
        as_tuple_t<std::make_integer_sequence<int, Rank>>());
  }

  template <typename State>
  State set_momentum(const State& q, const std::array<A, Rank>& rhou) const
      noexcept {
    using namespace variables;
    State result{q};
    std::array<A, Rank> old = get_momentum(q);
    auto e_kinetic_new = fub::transform_reduce(rhou, rhou, A(0));
    auto e_kinetic_old = fub::transform_reduce(old, old, A(0));
    result[energy] += 0.5 * (e_kinetic_new - e_kinetic_old) / q[density];
    for_each_tuple_element(
        [&](auto dir) { result[momentum<decltype(dir)::value>] = rhou[dir()]; },
        as_tuple_t<std::make_integer_sequence<int, Rank>>());
    return result;
  }

  template <typename Abi, typename State>
  std::enable_if_t<is_simd_abi_v<Abi>, State>
  set_momentum(const Abi& abi, const State& q,
               const nodeduce_t<std::array<simd<A, Abi>, Rank>>& rhou) const
      noexcept {
    using namespace variables;
    State result{q};
    std::array<simd<A, Abi>, Rank> old = get_momentum(abi, q);
    auto e_kinetic_new = fub::transform_reduce(rhou, rhou, simd<A, Abi>(0));
    auto e_kinetic_old = fub::transform_reduce(old, old, simd<A, Abi>(0));
    result[energy] += 0.5 * (e_kinetic_new - e_kinetic_old) / q[density];
    for_each_tuple_element(
        [&](auto dir) { result[momentum<decltype(dir)::value>] = rhou[dir()]; },
        as_tuple_t<std::make_integer_sequence<int, Rank>>());
    return result;
  }

  template <typename State>
  State set_velocity(const State& q, std::array<A, Rank> u) const noexcept {
    using variables::density;
    std::transform(u.begin(), u.end(), u.begin(),
                   [rho = q[density]](A u_i) { return rho * u_i; });
    return set_momentum(q, u);
  }

  template <typename Abi, typename State>
  std::enable_if_t<is_simd_abi_v<Abi>, State>
  set_velocity(const Abi& abi, const State& q,
               nodeduce_t<std::array<simd<A, Abi>, Rank>> u) const noexcept {
    using variables::density;
    std::transform(u.begin(), u.end(), u.begin(),
                   [rho = q[density]](simd<A, Abi> u_i) { return rho * u_i; });
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
  static std::array<A, species_size>
  get_mass_fractions(const State& q) noexcept {
    std::array<A, species_size> mass_fractions;
    fub::for_each_tuple_element(
        [&](auto species) {
          using variables::density;
          mass_fractions[as_index(species)] = q[species] / q[density];
        },
        tail_t<species_tuple>());
    mass_fractions[as_index(head_t<species_tuple>())] =
        std::max(A(0), 1 - std::accumulate(mass_fractions.begin() + 1,
                                           mass_fractions.end(), A(0)));
    return mass_fractions;
  }

  template <typename Abi, typename State>
  static std::array<simd<A, Abi>, species_size>
  get_mass_fractions(const Abi&, const State& q) noexcept {
    std::array<simd<A, Abi>, species_size> mass_fractions;
    fub::for_each_tuple_element(
        [&](auto species) {
          using variables::density;
          mass_fractions[as_index(species)] = q[species] / q[density];
        },
        tail_t<species_tuple>());
    mass_fractions[as_index(head_t<species_tuple>())] = fub::max(
        A(0), 1 - std::accumulate(mass_fractions.begin() + 1,
                                  mass_fractions.end(), simd<A, Abi>(0)));
    return mass_fractions;
  }
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // get mean specific enthalpy of formation from a state
  // {{{
  template <typename State>
  A get_specific_enthalpies_of_formation(const State& q) const noexcept {
    return q[variables::enthalpy];
  }
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // get mean specific heat capacity at constant pressure from a state
  // {{{
  A get_mean_specific_heat_capacity_at_constant_pressure(
      const complete_state& q) const noexcept {
    return q[variables::cp];
  }

  template <typename Abi>
  simd<A, Abi> get_mean_specific_heat_capacity_at_constant_pressure(
      Abi, const add_simd_t<complete_state, Abi>& q) const noexcept {
    return q[variables::cp];
  }

  A get_mean_specific_heat_capacity_at_constant_pressure(
      const conservative_state& q) const noexcept {
    using variables::temperature;
    auto cps =
        get_specific_heat_capacities_at_constant_pressure(q[temperature]);
    auto mass_fractions = get_mass_fractions(q);
    return transform_reduce(cps.begin(), cps.end(), mass_fractions.begin(),
                            0.0);
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
    auto mass_fractions = get_mass_fractions(Q);
    species_span<const A> Rs = get_specific_gas_constants();
    A R = fub::transform_reduce(Rs, mass_fractions, A(0));
    A h_ref = get_mean_reference_specific_enthalpy_of_formation(mass_fractions);
    A T_ref = get_reference_temperature();
    A U = Q[energy] / Q[density] + h_ref - R * T_ref;
    auto thermo = set_URY(U, Q[density], mass_fractions);
    return set_momentum(thermo, get_momentum(q));
  }

  template <typename Abi>
  add_simd_t<complete_state, Abi>
  get_complete_state(const Abi& abi,
                     const add_simd_t<conservative_state, Abi>& q) const {
    using namespace variables;
    add_simd_t<conservative_state, Abi> Q = set_momentum(abi, q, {0});
    auto mass_fractions = get_mass_fractions(abi, Q);
    species_span<const A> Rs = get_specific_gas_constants();
    simd<A, Abi> R = fub::transform_reduce(Rs, mass_fractions, simd<A, Abi>(0));
    simd<A, Abi> h_ref =
        get_mean_reference_specific_enthalpy_of_formation(abi, mass_fractions);
    simd<A, Abi> T_ref = get_reference_temperature();
    simd<A, Abi> U = Q[energy] / Q[density] + h_ref - R * T_ref;
    add_simd_t<complete_state, Abi> thermo =
        set_URY<Abi>(abi, U, Q[density], mass_fractions);
    return set_momentum(abi, thermo, get_momentum(abi, q));
  }
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // getter for variables
  // {{{
  template <typename Var> A get(Var, const complete_state& q) const noexcept {
    return q.template get<Var>();
  }

  template <int Dir>
  A get(Velocity<Dir>, const complete_state& q) const noexcept {
    using namespace variables;
    return q[momentum<Dir>] / q[density];
  }

  A get(head_t<species_tuple>, const complete_state& q) const noexcept {
    using namespace variables;
    A implicit_fraction = q[density] - foldl(tail_t<species_tuple>(), A(0),
                                             [&](A sum, auto species) {
                                               return sum + q[species];
                                             });
    return std::max(A(0), implicit_fraction);
  }

  template <typename Abi, typename Var>
  simd<A, Abi> get(Abi, Var, const add_simd_t<complete_state, Abi>& q) const
      noexcept {
    return q.template get<Var>();
  }

  template <typename Abi, int Dir>
  simd<A, Abi> get(Abi, Velocity<Dir>,
                   const add_simd_t<complete_state, Abi>& q) const noexcept {
    using namespace variables;
    return q[momentum<Dir>] / q[density];
  }

  template <typename Abi>
  simd<A, Abi> get(Abi, head_t<species_tuple>,
                   const add_simd_t<complete_state, Abi>& q) const noexcept {
    using namespace variables;
    simd<A, Abi> implicit_fraction =
        q[density] -
        foldl(tail_t<species_tuple>(), simd<A, Abi>(0),
              [&](simd<A, Abi> sum, auto species) { return sum + q[species]; });
    return fub::max(simd<A, Abi>(0), implicit_fraction);
  }
  // }}}
};

} // namespace euler
} // namespace fub

#endif // !FUB_EULER_IDEAL_GAS_HPP
