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

#ifndef FUB_EULER_PERFECT_GAS_HPP
#define FUB_EULER_PERFECT_GAS_HPP

#include "fub/euler/variables.hpp"
#include "fub/execution.hpp"
#include "fub/permutate_dimensions.hpp"
#include "fub/variable_list.hpp"
#include "fub/variable_view.hpp"

#include <cstring>

namespace fub {
inline namespace v1 {
namespace euler {

template <int Rank>
struct perfect_gas_complete
    : variable_list<Density, Momentum<Rank>, Energy, Pressure, SpeedOfSound> {
  using variable_list<Density, Momentum<Rank>, Energy, Pressure,
                      SpeedOfSound>::variable_list;
};

template <int Rank>
struct perfect_gas_conservative
    : variable_list<Density, Momentum<Rank>, Energy> {
  using variable_list<Density, Momentum<Rank>, Energy>::variable_list;
};

/// This equation class describes the single species perfect gas euler
/// equations.
///
/// The system is given by 3 equations
///
///    ρ_t    + (ρu)_x         = 0
///    (ρu)_t + (ρu^2 + p)_x   = 0
///    (ρE)_t + (u (ρE + p))_x = 0
///
/// and we consider the ideal gas equation of state
///
///    p = ρRT
template <int Rank, typename A = double> class perfect_gas {
public:
  using floating_point = A;

  using mdspan = basic_mdspan<floating_point, dynamic_extents_t<Rank>>;

  using const_mdspan =
      basic_mdspan<const floating_point, dynamic_extents_t<Rank>>;

  using complete = perfect_gas_complete<Rank>;

  using conservative = perfect_gas_conservative<Rank>;

  template <typename State> using view = basic_variable_view<State, mdspan>;

  template <typename State>
  using const_view = basic_variable_view<State, const_mdspan>;

  template <typename Accessor>
  using state_ref = variable_ref<complete, floating_point, Accessor>;

  template <typename Accessor>
  using state_cref = variable_ref<complete, const floating_point, Accessor>;

  template <typename Accessor>
  using cons_ref = variable_ref<conservative, floating_point, Accessor>;

  template <typename Accessor>
  using cons_cref = variable_ref<conservative, const floating_point, Accessor>;

  static constexpr int rank() noexcept { return Rank; }

  constexpr perfect_gas(floating_point gamma,
                        floating_point specific_gas_constant);

  constexpr floating_point gamma() const noexcept;

  constexpr floating_point specific_gas_constant() const noexcept;

  constexpr floating_point heat_capacity_at_constant_pressure() const noexcept;

  constexpr floating_point heat_capacity_at_constant_volume() const noexcept;

  template <typename Variable, typename Accessor>
  typename Accessor::value_type get(Variable var,
                                    const state_cref<Accessor>& state) const;

  template <std::ptrdiff_t Dim, typename Accessor>
  typename Accessor::value_type get(basic_tag<Velocity<Rank>, Dim>,
                                    const state_cref<Accessor>& state) const;

  template <typename AccessorCons, typename AccessorState>
  void reconstruct(const cons_cref<AccessorCons>& cons,
                   const state_ref<AccessorState>& state) const;

  template <typename AccessorLhs, typename AccessorRhs>
  void reconstruct(const state_ref<AccessorLhs>& cons,
                   const state_ref<AccessorRhs>& state) const;

  template <typename AccessorLhs, typename AccessorRhs>
  void reconstruct_prim(const state_ref<AccessorLhs>& prim,
                        const state_ref<AccessorRhs>& state) const;

  template <typename AccessorCons, typename AccessorState>
  void flux(const state_cref<AccessorCons>& state,
            const cons_ref<AccessorState>& flux) const;

  void permutate_dimensions(const_view<complete> from,
                            const std::array<int, 2>& permutation,
                            view<complete> to) {
    fub::permutate_dimensions(from, permutation, to);
    auto dyn_mom0 = dynamic_momentum<Rank>(permutation[0]);
    auto dyn_mom1 = dynamic_momentum<Rank>(permutation[1]);
    visit(
        [&](auto mom0, auto mom1) {
          std::swap_ranges(to[mom0].span().begin(), to[mom0].span().end(),
                           to[mom1].span().begin());
        },
        dyn_mom0, dyn_mom1);
  }

private:
  floating_point m_gamma;
  floating_point m_specific_gas_constant;
};

template <int Rank, typename A>
constexpr perfect_gas<Rank, A>::perfect_gas(A gamma, A specific_gas_constant)
    : m_gamma(gamma), m_specific_gas_constant(specific_gas_constant) {
  assert(m_gamma > A(1));
  assert(m_specific_gas_constant > A(0));
}

template <int Rank, typename A>
template <typename Variable, typename Accessor>
typename Accessor::value_type
perfect_gas<Rank, A>::get(Variable var,
                          const state_cref<Accessor>& state) const {
  return state[var];
}

template <int Rank, typename A>
template <std::ptrdiff_t Dim, typename Accessor>
typename Accessor::value_type
perfect_gas<Rank, A>::get(basic_tag<Velocity<Rank>, Dim>,
                          const state_cref<Accessor>& state) const {
  typename Accessor::value_type rho_u = state[momentum<Rank, Dim>];
  typename Accessor::value_type rho = state[density];
  assert(all_of(rho > 0));
  return rho_u / rho;
}

template <int Rank, typename A>
template <typename AccessorCons, typename AccessorState>
void perfect_gas<Rank, A>::reconstruct(
    const cons_cref<AccessorCons>& cons,
    const state_ref<AccessorState>& state) const {
  static_assert(std::is_same<typename AccessorState::value_type,
                             typename AccessorCons::value_type>{},
                "Inconsistency: Accessor have different value type. This is "
                "probably an error.");
  using float_t = typename AccessorState::value_type;
  const float_t gamma_minus_1 = m_gamma - 1.f;
  for_each(cons.get_variable_list(),
           [&](auto variable) { state[variable] = cons[variable]; });
  const float_t rho = cons[density];
  const float_t rho_e = cons[energy];
  const float_t e_kin = [&] {
    float_t e_kin = A{0};
    for_each(Momentum<Rank>(), [&](auto mom) {
      const float_t rho_u = state[mom];
      e_kin += rho_u * rho_u;
    });
    e_kin *= float_t(0.5) / rho;
    return rho * e_kin;
  }();
  const float_t e_int = rho_e - e_kin;
  state[pressure] = gamma_minus_1 * e_int;
  state[speed_of_sound] = std::sqrt(m_gamma * state[pressure] / state[density]);
}

template <int Rank, typename A>
template <typename AccessorLhs, typename AccessorRhs>
void perfect_gas<Rank, A>::reconstruct(
    const state_ref<AccessorLhs>& cons,
    const state_ref<AccessorRhs>& state) const {
  static_assert(std::is_same<typename AccessorLhs::value_type,
                             typename AccessorRhs::value_type>{},
                "Inconsistency: Accessor have different value type. This is "
                "probably an error.");
  using float_t = typename AccessorLhs::value_type;
  const float_t gamma_minus_1 = m_gamma - 1.f;
  for_each(cons.get_variable_list(),
           [&](auto variable) { state[variable] = cons[variable]; });
  const float_t rho = cons[density];
  const float_t rho_e = cons[energy];
  const float_t e_kin = [&] {
    float_t e_kin = A{0};
    for_each(Momentum<Rank>(), [&](auto mom) {
      const float_t rho_u = state[mom];
      e_kin += rho_u * rho_u;
    });
    e_kin *= float_t(0.5) / rho;
    return rho * e_kin;
  }();
  const float_t e_int = rho_e - e_kin;
  state[pressure] = gamma_minus_1 * e_int;
  state[speed_of_sound] = std::sqrt(m_gamma * state[pressure] / state[density]);
}

template <int Rank, typename A>
template <typename AccessorLhs, typename AccessorRhs>
void perfect_gas<Rank, A>::reconstruct_prim(
    const state_ref<AccessorLhs>& prim,
    const state_ref<AccessorRhs>& state) const {
  static_assert(std::is_same<typename AccessorLhs::value_type,
                             typename AccessorRhs::value_type>{},
                "Inconsistency: Accessor have different value type. This is "
                "probably an error.");
  using float_t = typename AccessorLhs::value_type;
  const float_t gamma_minus_1 = m_gamma - 1.f;
  for_each(prim.get_variable_list(),
           [&](auto variable) { state[variable] = prim[variable]; });
  const float_t rho = prim[density];
  const float_t p = prim[pressure];
  const float_t e_kin = [&] {
    float_t e_kin = A{0};
    for_each(Momentum<Rank>(), [&](auto mom) {
      const float_t rho_u = state[mom];
      e_kin += rho_u * rho_u;
    });
    e_kin *= float_t(0.5) / rho;
    return rho * e_kin;
  }();
  const float_t e_int = p / gamma_minus_1;
  state[energy] = e_kin + e_int;
  state[speed_of_sound] = std::sqrt(m_gamma * p / rho);
}

template <int Rank, typename A>
template <typename AccessorState, typename AccessorCons>
void perfect_gas<Rank, A>::flux(const state_cref<AccessorState>& state,
                                const cons_ref<AccessorCons>& flux) const {
  static_assert(std::is_same<typename AccessorState::value_type,
                             typename AccessorCons::value_type>{},
                "Inconsistency: Accessor have different value type. This is "
                "probably an error.");
  using float_t = typename AccessorState::value_type;
  const float_t rho = state[density];
  assert(all_of(rho > 0));
  const float_t rho_u = state[momentum<Rank, 0>];
  const float_t p = state[pressure];
  const float_t rho_e = state[energy];
  for_each(Momentum<Rank>(), [&](auto momentum) {
    const float_t rho_v = state[momentum];
    flux[momentum] = rho_u / rho * rho_v;
  });
  flux[density] = rho_u;
  flux[momentum<Rank, 0>] += p;
  flux[energy] = rho_u / rho * (rho_e + p);
}

template <int Rank, typename A>
constexpr A perfect_gas<Rank, A>::heat_capacity_at_constant_volume() const
    noexcept {
  return m_specific_gas_constant / (m_gamma - 1);
}

} // namespace euler
} // namespace v1
} // namespace fub

#endif