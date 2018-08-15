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
#include "fub/span.hpp"
#include "fub/tuple.hpp"

#include <numeric>
#include <tuple>

namespace fub {
namespace euler {

/// This problem describes the multiple species ideal gas euler equations.
///
/// The system is given by s + 3 equations for s + 1 species
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
template <typename Mechanism, int Rank = 1> class ideal_gas {
public:
  /// \name Types
  using floating_point = typename Mechanism::floating_point;

  /// This defines all variables on a data patch for this problem.
  struct complete_state_variables
      : variable_list<Density, Energy, Pressure, Temperature, SpeedOfSound,
                      HeatCapacityAtConstantP, Enthalpy, Momentum<Rank>,
                      Species> {
    using variable_list<Density, Energy, Pressure, Temperature, SpeedOfSound,
                        HeatCapacityAtConstantP, Enthalpy, Velocity<Rank>,
                        Species>::variable_list;
  };

  using complete_state = quantities<complete_state_variables>;

  /// This defines all variables which shall be conserved by the finite volume
  /// methods.
  struct conservative_state_variables
      : variable_list<Density, Momentum<Rank>, Energy, Species> {};

  using conservative_state = quantities<conservative_state_variables>;

  /// \name Mechanism Functions

  /// Returns an array of molar masses for each species.
  span<const floating_point> get_molar_masses() const noexcept {
    return m_mechanism.get_molar_masses();
  }

  /// Computes the production rates for each species and stores them in dYdt.
  ///
  /// \param[out] dYdt  The production rate results.
  ///
  /// \param[in] Y  The current mass fractions for each species.
  ///
  /// \param[in] T  A temperature value in kelvin [K].
  ///
  /// \param[in] P  A pressure value in Pascal [Pa].
  void get_production_rates(span<floating_point> dYdt,
                            span<const floating_point> Y, floating_point T,
                            floating_point P) const noexcept {
    m_mechanism.get_production_rates(dYdt, Y, T, P);
  }

  void get_specific_enthalpies_of_formation(span<floating_point> enthalpies,
                                            floating_point temperature) const
      noexcept {
    m_mechanism.get_specific_enthalpies_of_formation(enthalpies, temperature);
  }

  void get_specific_heat_capacities_at_constant_pressure(
      span<floating_point> heat_capacities, floating_point temperature) const
      noexcept {
    m_mechanism.get_specific_heat_capacities_at_constant_pressure(
        heat_capacities, temperature);
  }

private:
  /// Returns an array which contains all specific gas constants.
  auto store_specific_gas_constants() const noexcept {
    // [J / (kmol K)]
    const floating_point universal_gas_constant = get_universal_gas_constant();
    // [kg / kmol]
    span<const floating_point> molar_masses = get_molar_masses();

    std::array<floating_point, species_size> constants;
    std::transform(molar_masses.begin(), molar_masses.end(), constants.begin(),
                   [=](floating_point molar_mass) {
                     // [J / (kg K)]
                     return universal_gas_constant / molar_mass;
                   });
    return constants;
  }

public:
  /// Returns a view to a static array containing specific gas constants
  /// for each species.
  span<const floating_point> get_specific_gas_constants() const noexcept {
    static const auto specific_gas_constants = store_specific_gas_constants();
    return specific_gas_constants;
  }

  constexpr floating_point get_reference_temperature() const noexcept {
    return floating_point(300);
  }

private:
  auto store_reference_specific_enthalpies_of_formation() const noexcept {
    const floating_point T_ref = get_reference_temperature();

    std::array<floating_point, species_size> constants;
    get_specific_enthalpies_of_formation(constants, T_ref);
    return constants;
  }

public:
  /// Returns the specific enthalpy at reference temeprature for the
  /// specified mass fractions.
  floating_point get_mean_reference_specific_enthalpy_of_formation(
      span<const floating_point> mass_fractions) const noexcept {
    static const auto h0 = store_reference_specific_enthalpies_of_formation();
    return fub::transform_reduce(h0, mass_fractions, floating_point(0));
  }

  /// Returns the universal gas constant which this mechanism assumes.
  floating_point get_universal_gas_constant() const noexcept {
    return floating_point(8314.45);
  }

  /// Computes the internal energy dependent on temperature, enthalpy
  /// and mass fractions.
  ///
  /// We are using the formula
  ///
  ///     U = \int c_v(T) dT = \int (c_p(T) - R) dT = H - RT.
  ///
  /// Where H is the mean enthalpy, R the specific gas constant and `T`
  /// temperature.
  floating_point get_specific_internal_energy(
      floating_point temperature, floating_point mean_specific_gas_constant,
      span<const floating_point> enthalpies,
      span<const floating_point> mass_fractions) const noexcept {
    // [J / kg]
    floating_point mean_H =
        fub::transform_reduce(enthalpies, mass_fractions, floating_point(0));
    return mean_H - mean_specific_gas_constant * temperature;
  }

private:
  // Returns the temperature which matches a specified intenral
  // energy.
  //
  // This function uses a Newton-like iteration to find a temperature such
  // that
  //
  //     get_internal_energy(T, Y_0) == U_0
  //
  // for specified Y_0 and U_0.
  floating_point find_temperature(span<floating_point> h,
                                  span<floating_point> cps, floating_point U_0,
                                  span<const floating_point> Y_0,
                                  floating_point T_0 = 300,
                                  std::ptrdiff_t tolerance = 100000,
                                  std::ptrdiff_t max_iterations = 1000) const {
    const span<const floating_point> Rs = get_specific_gas_constants();
    const floating_point R = fub::transform_reduce(Rs, Y_0, A(0));
    auto get_U = [&](floating_point T) {
      get_specific_enthalpies_of_formation(h, T);
      get_specific_heat_capacities_at_constant_pressure(cps, T);
      return get_specific_internal_energy(T, R, h, Y_0);
    };
    floating_point T = T_0;
    floating_point U = get_U(T);
    std::ptrdiff_t iterations = 0;
    while (!almost_equal(U_0, U, tolerance) && iterations < max_iterations) {
      const floating_point cp = fub::transform_reduce(cps, Y_0, A(0));
      const floating_point cv = cp - R;
      const floating_point dU = U_0 - U;
      const floating_point dT = fub::clamp(dU / cv, A(-100), A(100));
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

public:
  /// Returns a valid and complete set of variables computed from
  /// conservative variables alone.
  ///
  /// This function will be used to retain a complete state after integrating
  /// the conservative variables in time
  complete_state set_URY(floating_point internal_energy, floating_point rho,
                         span<const floating_point> mass_fractions) const {
    using namespace variables;
    complete_state q;
    q[density] = rho;
    std::vector<floating_point> hs(species().size());
    std::vector<floating_point> cps(species().size());
    q[temperature] = find_temperature(hs, cps, internal_energy, mass_fractions);
    q[cp] = fub::transform_reduce(cps, mass_fractions, floating_point(0));
    q[enthalpy] = fub::transform_reduce(hs, mass_fractions, floating_point(0));
    span<const floating_point> specific_gas_constants =
        get_specific_gas_constants();
    const floating_point R = fub::transform_reduce(
        specific_gas_constants, mass_fractions, floating_point(0));
    q[pressure] = rho * R * q[temperature];
    const floating_point gamma = q[cp] / (q[cp] - R);
    q[speed_of_sound] = std::sqrt(gamma * R * q[temperature]);
    hana::for_each(species(), [&](auto species) {
      q[species] = rho * mass_fractions[species.index()];
    });
    const floating_point h_ref =
        get_mean_reference_specific_enthalpy_of_formation(mass_fractions);
    const floating_point T_ref = get_reference_temperature();
    q[energy] = rho * (internal_energy + R * T_ref - h_ref);
    return q;
  }

  void
  get_mass_from_mole_fractions(span<floating_point> fractions,
                               span<const floating_point> mole_fractions) const
      noexcept {
    span<const floating_point> molar_masses = get_molar_masses();
    std::transform(mole_fractions.begin(), mole_fractions.end(),
                   molar_masses.begin(), fractions.begin(),
                   std::multiplies<>{});
    std::transform(fractions.begin(), fractions.end(), fractions.begin(),
                   [](double x) { return std::max(floating_point(0), x); });
    floating_point sum =
        std::accumulate(fractions.begin(), fractions.end(), floating_point(0));
    std::transform(fractions.begin(), fractions.end(), fractions.begin(),
                   [sum](floating_point Y) {
                     return std::min(floating_point(1), Y / sum);
                   });
    return fractions;
  }

  complete_state set_TPY(span<floating_point> hs, span<floating_point> cps,
                         floating_point T, floating_point p,
                         span<const floating_point> mass_fractions) const
      noexcept {
    using namespace variables;
    complete_state q;
    q[temperature] = T;
    q[pressure] = p;
    span<const floating_point> Rs = get_specific_gas_constants();
    const floating_point R =
        fub::transform_reduce(Rs, mass_fractions, floating_point(0));
    q[density] = p / (R * T);
    {
      // Set total energy
      get_specific_enthalpies_of_formation(hs, T);
      q[enthalpy] = fub::transform_reduce(hs, mass_fractions, A(0));
      const floating_point U =
          get_specific_internal_energy(T, R, hs, mass_fractions);
      const floating_point T_ref = get_reference_temperature();
      const floating_point h_ref =
          get_mean_reference_specific_enthalpy_of_formation(mass_fractions);
      const floating_point E_ref = R * T_ref - h_ref;
      q[energy] = q[density] * (U + E_ref);
    }
    get_specific_heat_capacities_at_constant_pressure(cps, T);
    q[cp] = fub::transform_reduce(cps, mass_fractions, A(0));
    const floating_point gamma = q[cp] / (q[cp] - R);
    assert(T >= 0);
    assert(gamma >= 0);
    assert(R >= 0);
    q[speed_of_sound] = std::sqrt(gamma * R * T);
    floating_point species_sum = 0;
    hana::for_each(species(), [&](auto species) {
      q[species] = rho * mass_fractions[species.index()];
    });
    return q;
  }

  complete_state set_TPY(floating_point T, floating_point p,
                         span<const floating_point> mass_fractions) const
      noexcept {
    std::vector<floating_point> hs(species().size());
    std::vector<floating_point> cps(species().size());
    return set_TPY(hs, cps, T, p, mass_fractions);
  }

  complete_state set_TPX(span<floating_point> hs, span<floating_point> cps,
                         floating_point temperature, floating_point pressure,
                         span<const floating_point> mole_fractions) const
      noexcept {
    std::vector<floating_point> mass_fractions(species().size());
    get_mass_from_mole_fractions(mass_fractions, mole_fractions);
    return set_TPY(hs, cps, temperature, pressure, mass_fractions);
  }

  complete_state set_TPX(floating_point temperature, floating_point pressure,
                         span<const floating_point> mole_fractions) const
      noexcept {
    std::vector<floating_point> mass_fractions(species().size());
    get_mass_from_mole_fractions(mass_fractions, mole_fractions);
    return set_TPY(temperature, pressure, mass_fractions);
  }

  template <typename State>
  std::array<floating_point, Rank> get_momentum(const State& q) const noexcept {
    std::array<floating_point, Rank> array;
    for_each(Momentum<Rank>(),
             [&](auto momentum) { array[momentum.index()] = q[momentum]; });
    return q;
  }

  template <typename State>
  State set_momentum(const State& q,
                     span<const floating_point, Rank> rhou) const noexcept {
    using namespace variables;
    State result{q};
    std::array<floating_point, Rank> old = get_momentum(q);
    auto e_kinetic_new = fub::transform_reduce(rhou, rhou, floating_point(0));
    auto e_kinetic_old = fub::transform_reduce(old, old, floating_point(0));
    result[energy] += 0.5 * (e_kinetic_new - e_kinetic_old) / q[density];
    for_each(Momentum<Rank>(),
             [&](auto momentum) { result[momentum] = rhou[momentum.index()]; });
    return result;
  }

  template <typename State>
  State set_velocity(const State& q, span<const floating_point, Rank> u) const
      noexcept {
    std::transform(
        u.begin(), u.end(), u.begin(),
        [rho = q[density]](floating_point u_i) { return rho * u_i; });
    return set_momentum(q, u);
  }

  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // get flux
  // {{{

  /// Returns and evaluates the euler flux function for conservative
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
    for_each_tuple_element(
        [&](auto fq) { f[fq] = u * q[unflux(fq)]; },
        fub::apply([](auto... vs) { return flatten_variables(vs...); },
                   get_variables(f)));
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
    for_each_tuple_element(
        [&](auto fq) { f[fq] = u * q[unflux(fq)]; },
        fub::apply([](auto... vs) { return flatten_variables(vs...); },
                   get_variables(f)));
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
    mass_fractions[as_index(head_t<species_tuple>())] =
        fub::max(simd<A, Abi>(A(0)),
                 1 - std::accumulate(mass_fractions.begin() + 1,
                                     mass_fractions.end(), simd<A, Abi>(0)));
    return mass_fractions;
  }
  // }}}

  //////////////////////////////////////////////////////////////////////////////
  // get mean specific heat capacity at constant pressure from a state
  // {{{
  floating_point get_mean_specific_heat_capacity_at_constant_pressure(
      const complete_state& q) const noexcept {
    return q[cp];
  }

  floating_point get_mean_specific_heat_capacity_at_constant_pressure(
      const conservative_state& q) const noexcept {
    auto cps =
        get_specific_heat_capacities_at_constant_pressure(q[temperature]);
    auto mass_fractions = get_mass_fractions(q);
    return fub::transform_reduce(cps.begin(), cps.end(), mass_fractions.begin(),
                                 0.0);
  }
  // }}}

  /// Retains a complete state from a conservative one.
  ///
  /// This function maybe used after a state update of conservative variables.
  complete_state get_complete_state(const conservative_state& q) const {
    using namespace variables;
    conservative_state Q = set_momentum(q, {0});
    auto mass_fractions = get_mass_fractions(Q);
    span<const A> Rs = get_specific_gas_constants();
    A R = fub::transform_reduce(Rs, mass_fractions, A(0));
    A h_ref = get_mean_reference_specific_enthalpy_of_formation(mass_fractions);
    A T_ref = get_reference_temperature();
    A U = Q[energy] / Q[density] + h_ref - R * T_ref;
    auto thermo = set_URY(U, Q[density], mass_fractions);
    return set_momentum(thermo, get_momentum(q));
  }

  template <typename Var> A get(Var, const complete_state& q) const noexcept {
    return q.template get<Var>();
  }

  template <typename Archive> void serialize(Archive& archive, unsigned) {
    archive& m_mechanism;
  }

private:
#if __has_cpp_attribute(no_unique_address)
  [[no_unique_address]] Mechanism m_mechanism;
#else
  Mechanism m_mechanism;
#endif
};

} // namespace euler
} // namespace fub

#endif // !FUB_EULER_IDEAL_GAS_HPP
