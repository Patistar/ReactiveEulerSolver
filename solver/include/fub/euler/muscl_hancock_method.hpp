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

#ifndef FUB_EULER_MUSCL_HANCOCK_METHOD_HPP
#define FUB_EULER_MUSCL_HANCOCK_METHOD_HPP

#include "fub/equation.hpp"
#include "fub/euler/ideal_gas.hpp"
#include "fub/godunov_method.hpp"
#include "fub/patch_view.hpp"
#include "fub/time_integrator/forward_euler.hpp"
#include "fub/tuple.hpp"

#include <chrono>
#include <tuple>

namespace fub {
namespace euler {

/// @brief This method is a MUSCL-type flux method which deploys a limited state
/// reconstruction on a halftime-level.
///
/// State reconstruction happens for the primitve variables `T` (temperature),
/// `u` (velocity), `p` (pressure) and `{Y_s, ...}` (mass fractions for all
/// species `s`).
template <typename RiemannSolver> struct muscl_hancock_method {
  godunov_method<RiemannSolver> m_godunov_method{};

  /// @brief The default constructor depends on godunov_method's default
  /// constructor
  muscl_hancock_method() = default;

  /// @brief Initialises the internal godunov_method with the specified riemann
  /// problem solver.
  explicit muscl_hancock_method(const RiemannSolver& rps)
      : m_godunov_method(rps) {}

  /// @brief Directly intialises this method by the godunov method.
  explicit muscl_hancock_method(const godunov_method<RiemannSolver>& method)
      : m_godunov_method(method) {}

  template <int Rank> struct velocity_states {
    template <typename> struct impl;
    template <int... Is> struct impl<std::integer_sequence<int, Is...>> {
      using type = quantities<Velocity<Is>...>;
    };
    using type = typename impl<std::make_integer_sequence<int, Rank>>::type;
  };
  template <int Rank>
  using velocity_states_t = typename velocity_states<Rank>::type;

  /// @brief Returns the primtive state vector for a specified Mechanism
  ///
  /// The state vector contains velocity, pressure and temperature plus all
  /// relative species mass fractions, including the first species.
  template <typename Mechanism, int Rank>
  using primitive_state_t =
      concat_t<velocity_states_t<Rank>, quantities<Pressure, Temperature>,
               as_quantities_t<species_tuple_t<ideal_gas<Mechanism, Rank>>>>;

  /// @brief Returns all primtive variables as a tuple
  template <typename Mechanism, int Rank>
  using primitive_tuple_t = as_tuple_t<primitive_state_t<Mechanism, Rank>>;

  /// @brief Returns a simd-enabled primtive state vector type for the specified
  /// simd abi.
  template <typename Mechanism, int Rank, typename Abi>
  using primitive_simd_t = add_simd_t<primitive_state_t<Mechanism, Rank>, Abi>;

  /// @brief We use this function as a building block to compute limited slopes
  template <typename T, typename Abi>
  static simd<T, Abi> minmod(const simd<T, Abi>& a,
                             const nodeduce_t<simd<T, Abi>>& b) noexcept {
    simd<T, Abi> m = a;
    where(abs(b) < abs(a), m) = b;
    where(a * b <= 0, m) = 0;
    return m;
  }

  /// @brief We deploy Monotonized central-difference (MC) Limiter
  template <typename T, typename Abi>
  static simd<T, Abi>
  compute_limited_slope(const simd<T, Abi>& left,
                        const nodeduce_t<simd<T, Abi>>& right) noexcept {
    simd<T, Abi> r = right / left;
    simd<T, Abi> van_leer_limiter  = (4 * r) / ((r + 1) * (r + 1));
    where(left * right <= 0, van_leer_limiter) = 0;
    return van_leer_limiter * left;
  }

  /// @brief Returns the slope vector for a specified stencil
  template <typename Mechanism, int Rank, typename Abi>
  static auto
  compute_slopes(const ideal_gas<Mechanism, Rank>&, const Abi&,
                 const primitive_simd_t<Mechanism, Rank, Abi>& left,
                 const primitive_simd_t<Mechanism, Rank, Abi>& mid,
                 const primitive_simd_t<Mechanism, Rank, Abi>& right) noexcept {
    primitive_simd_t<Mechanism, Rank, Abi> slope;
    for_each_tuple_element(
        [&](auto var) {
          const auto dQ_left = mid[var] - left[var];
          const auto dQ_right = right[var] - mid[var];
          const auto limited_slope = compute_limited_slope(dQ_left, dQ_right);
          slope[var] = limited_slope;
        },
        get_variables(slope));
    return slope;
  }

  /// @brief Returns a state vector containing all primitive variables.
  ///
  /// If neccessary this function will compute the variables if they are not
  /// directly available in the specified state.
  ///
  /// Note that the equation knows how to compute them.
  template <typename Mechanism, int Rank, typename Abi, typename State>
  static primitive_simd_t<Mechanism, Rank, Abi>
  as_primitive_state(const ideal_gas<Mechanism, Rank>& eq, const Abi& abi,
                     const State& q) noexcept {
    using namespace variables;
    primitive_simd_t<Mechanism, Rank, Abi> w;
    w[pressure] = eq.get(abi, pressure, q);
    w[temperature] = eq.get(abi, temperature, q);
    for_each_tuple_element(
        [&](auto velocity) { w[velocity] = eq.get(abi, velocity, q); },
        as_tuple_t<velocity_states_t<Rank>>());
    auto rho = eq.get(abi, density, q);
    for_each_tuple_element(
        [&](auto species) { w[species] = eq.get(abi, species, q) / rho; },
        species_tuple(eq));
    return w;
  }

  /// @brief We use this pair-like type to give the states some names.
  template <typename Mechanism, int Rank, typename Abi> struct evolved_states {
    using state_type =
        add_simd_t<complete_state_t<ideal_gas<Mechanism, Rank>>, Abi>;
    state_type left;
    state_type right;
  };

  /// @brief Returns two complete states located at left and right face
  /// boundaries for the halftime-level. This function expects a cell centered
  /// base state with specified limited slope `dwdx`.
  template <typename Mechanism, int Rank, typename Abi>
  evolved_states<Mechanism, Rank, Abi>
  evolve(const ideal_gas<Mechanism, Rank>& equation, const Abi& abi,
         const primitive_simd_t<Mechanism, Rank, Abi>& w,
         const nodeduce_t<simd<double, Abi>>& rho,
         const nodeduce_t<simd<double, Abi>>& cp,
         const primitive_simd_t<Mechanism, Rank, Abi>& dwdx,
         double lambda) const noexcept {
    using namespace variables;
    using A = double;
    using S = simd<A, Abi>;
    static constexpr index size = ideal_gas<Mechanism, Rank>::species_size;
    std::array<S, size> mass_fractions;
    fub::for_each_tuple_element(
        [&](auto species) { mass_fractions[as_index(species)] = w[species]; },
        species_tuple(equation));
    span<const A, size> Rs = equation.get_specific_gas_constants();
    S R = fub::transform_reduce(Rs, mass_fractions, S(0));
    S R_hat = equation.get_universal_gas_constant();
    S cv = cp - R;

    // Compute Change in Species Concentrations, dY / dt
    std::array<S, size> dYdt;
    fub::for_each_tuple_element(
        [&](auto species) { dYdt[species] = -w[velocity<0>] * dwdx[species]; });

    // Compute, dR / dt
    auto molar_masses = equation.get_molar_masses();
    S dRdt = R_hat * fub::transform_reduce(dYdt, molar_masses, S(0),
                                           std::plus<>(), std::divides<>());

    // Compute spatial density gradient, drho / dx
    using namespace variables;
    S dYdx_mass_sum = 0;
    fub::for_each_tuple_element(
        [&](auto species) {
          dYdx_mass_sum += dwdx[species] / molar_masses[as_index(species)];
        },
        species_tuple(equation));
    S drhodx = dwdx[pressure] / (R * w[temperature]) -
               w[pressure] / (R * w[temperature] * w[temperature]) *
                   dwdx[temperature] -
               (rho / R) * R_hat * dYdx_mass_sum;

    // Compute velocity derivate, du / dt
    S dudt = -w[velocity<0>] * dwdx[velocity<0>] - dwdx[pressure] / rho;

    // Compute temperature derivative, dT / dt
    S dTdt = -(w[pressure] * dwdx[velocity<0>] +
               rho * w[velocity<0>] * cv * dwdx[temperature] +
               0.5 * w[velocity<0>] * dwdx[pressure]) /
             (rho * cv);

    // Compute pressure derivate, dp / dt
    // dpdt = -(drhodx .* u + rho .* dudx) .* R .* T + rho .* R .* dTdt + rho .*
    // dRdt .* T;
    S dpdt = -(drhodx * w[velocity<0>] + rho * dwdx[velocity<0>]) * R *
                 w[temperature] +
             rho * R * dTdt + rho * dRdt * w[temperature];

    // Construct States
    using species_t = species_tuple_t<ideal_gas<Mechanism, Rank>>;
    static constexpr std::size_t array_size = std::tuple_size<species_t>::value;
    std::array<S, array_size> Y;
    evolved_states<Mechanism, Rank, Abi> states;

    // Construct left state
    S p = w[pressure] + 0.5 * (lambda * dpdt - dwdx[pressure]);
    std::array<S, Rank> u;
    u[0] = w[velocity<0>] + 0.5 * (lambda * dudt - dwdx[velocity<0>]);
    for_each_tuple_element(
        [&](auto dim) {
          static constexpr int Dim = decltype(dim)::value;
          auto dvdt = -w[velocity<Dim>] * dwdx[velocity<Dim>];
          u[Dim] =
              w[velocity<Dim>] + 0.5 * (lambda * dvdt - dwdx[velocity<Dim>]);
        },
        tail_t<as_tuple_t<std::make_integer_sequence<int, Rank>>>());
    S T = w[temperature] + 0.5 * (lambda * dTdt - dwdx[temperature]);
    for_each_tuple_element(
        [&](auto species) {
          Y[as_index(species)] =
              w[species] +
              0.5 * (lambda * dYdt[as_index(species)] - dwdx[species]);
        },
        species_tuple(equation));
    states.left = equation.set_velocity(abi, equation.set_TPY(abi, T, p, Y), u);

    // Construct right state
    p = w[pressure] + 0.5 * (lambda * dpdt + dwdx[pressure]);
    u[0] = w[velocity<0>] + 0.5 * (lambda * dudt + dwdx[velocity<0>]);
    for_each_tuple_element(
        [&](auto dim) {
          static constexpr int Dim = decltype(dim)::value;
          auto dvdt = -w[velocity<Dim>] * dwdx[velocity<Dim>];
          u[Dim] =
              w[velocity<Dim>] + 0.5 * (lambda * dvdt + dwdx[velocity<Dim>]);
        },
        tail_t<as_tuple_t<std::make_integer_sequence<int, Rank>>>());
    T = w[temperature] + 0.5 * (lambda * dTdt + dwdx[temperature]);
    for_each_tuple_element(
        [&](auto species) {
          Y[as_index(species)] =
              w[species] +
              0.5 * (lambda * dYdt[as_index(species)] + dwdx[species]);
        },
        species_tuple(equation));
    states.right =
        equation.set_velocity(abi, equation.set_TPY(abi, T, p, Y), u);

    return states;
  }

  template <typename Mechanism, int Rank, index Size>
  struct reconstructed_states {
    reconstructed_states(const extents<Size>& e)
        : on_left_face(e), on_right_face(e) {}

    static_assert(Size == dyn || Size > 0, "Invalid Extents Size.");
    using patch_type =
        patch<complete_state_t<ideal_gas<Mechanism, Rank>>, extents<Size>>;

    patch_type on_left_face;
    patch_type on_right_face;
  };

  template <typename Mechanism, int Rank, typename Abi, typename StateSpan,
            typename Coordinates>
  std::enable_if_t<
      is_simd_abi_v<Abi>,
      reconstructed_states<Mechanism, Rank, view_extents_t<StateSpan>().static_get(0)>>
  reconstruct_states(const ideal_gas<Mechanism, Rank>& eq, const Abi& abi,
                     const StateSpan& qL, const nodeduce_t<StateSpan>& qM,
                     const nodeduce_t<StateSpan>& qR,
                     std::chrono::duration<double> dt,
                     const Coordinates& coordinates) const noexcept {
    static constexpr index Extent = view_extents_t<StateSpan>().static_get(0);
    reconstructed_states<Mechanism, Rank, Extent> rec(qM.extents());
    double lambda = dt.count() / coordinates.dx()[0];
    fub::for_each_simd(
        abi,
        [&](auto abi, auto&& left, auto&& right, auto qL, auto qM, auto qR) {
          using prim_t = primitive_simd_t<Mechanism, Rank, decltype(abi)>;
          prim_t wL = as_primitive_state(eq, abi, qL);
          prim_t wM = as_primitive_state(eq, abi, qM);
          prim_t wR = as_primitive_state(eq, abi, qR);
          auto slope = compute_slopes(eq, abi, wL, wM, wR);
          auto cpM =
              eq.get_mean_specific_heat_capacity_at_constant_pressure(abi, qM);
          auto rho = eq.get(abi, variables::density, qM);
          auto states = evolve(eq, abi, wM, rho, cpM, slope, lambda);
          left = states.left;
          right = states.right;
        },
        make_view(rec.on_left_face), make_view(rec.on_right_face), qL, qM, qR);
    return rec;
  }

  template <typename Abi, typename Equation, typename Coordinates, typename F,
            typename L, typename M, typename R>
  void compute_numeric_fluxes(const Abi& abi, const F& flux, const L& left,
                              const M& middle, const R& right,
                              std::chrono::duration<double> dt,
                              const Coordinates& coordinates,
                              const Equation& equation) const {
    static_assert(is_view<L>::value, "The left data is not a view type.");
    static_assert(is_view<M>::value, "The middle data is not a view type.");
    static_assert(is_view<R>::value, "The right data is not a view type.");

    fub::for_each_row(
        [&](auto flux, auto left, auto middle, auto right) {
          // Gather all states into one continuous array
          const auto row = join(rtake<2>(left), middle, take<2>(right));
          const auto rowv = make_view(row);
          auto qL = rdrop<2>(rowv);
          auto qM = rdrop<1>(drop<1>(rowv));
          auto qR = drop<2>(rowv);

          // Reconstruct states on half-time level
          const auto reconstructed = reconstruct_states(
              equation, abi, qL, qM, qR, dt, coordinates);

          // Compute numeric fluxes with these reconstruced states.
          m_godunov_method.compute_numeric_fluxes(
              abi, flux, rdrop<1>(make_view(reconstructed.on_right_face)),
              drop<1>(make_view(reconstructed.on_left_face)), dt, coordinates,
              equation);
        },
        flux, left, middle, right);
  }

  template <typename Equation, typename Coordinates, typename F, typename L,
            typename M, typename R>
  void compute_numeric_fluxes(const F& flux, const L& left, const M& middle,
                              const R& right, std::chrono::duration<double> dt,
                              const Coordinates& coordinates,
                              const Equation& equation) const {
    compute_numeric_fluxes(simd_abi::native<double>(), flux, left, middle,
                           right, dt, coordinates, equation);
  }

  template <typename Equation, typename L, typename M, typename R,
            typename Coordinates>
  std::enable_if_t<conjunction<is_view<L>, is_view<R>, is_view<M>>::value,
                   std::chrono::duration<double>>
  get_stable_time_step(const Equation& equation, const L& left, const M& middle,
                       const R& right, const Coordinates& coordinates) const {
    return m_godunov_method.get_stable_time_step(equation, left, middle, right,
                                                 coordinates);
  }

  template <typename Archive>
  void serialize(Archive& ar, unsigned) {
      ar & m_godunov_method;
  }
};

} // namespace euler
} // namespace fub

#endif // !FUB_EULER_MUSCL_HANCOCK_METHOD_HPP
