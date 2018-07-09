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

#ifndef FUB_EULER_KINETICSOURCETERM_HPP
#define FUB_EULER_KINETICSOURCETERM_HPP

#include "fub/euler/ideal_gas.hpp"
#include "fub/ode_solver/radau.hpp"
#include "fub/patch_view.hpp"
#include "fub/serial/grid.hpp"
#include "fub/variables.hpp"

#include <array>
#include <chrono>
#include <numeric>

namespace fub {
namespace euler {
namespace kinetic_source_term_detail {
template <typename Grid, typename Solver> struct advance_patch;
template <typename Grid, typename Solver> struct advance_patch_action;
template <typename Grid, typename Solver> struct get_dt;
template <typename Grid, typename Solver> struct get_dt_action;
} // namespace kinetic_source_term_detail

template <typename Grid, typename OdeSolver = ode_solver::Radau>
struct kinetic_source_term {
  using traits = grid_traits<Grid>;
  using equation_type = typename traits::equation_type;
  using partition_type = typename traits::partition_type;
  using const_patch_view = typename traits::const_patch_view_type;
  using patch_type = typename traits::patch_type;
  using node_type = typename traits::node_type;

  template <typename T> using future = typename traits::template future<T>;
  using location_type = typename traits::location_type;

  using complete_state = typename equation_type::complete_state;

  using duration = std::chrono::duration<double>;

  static constexpr int Size =
      std::tuple_size<typename equation_type::species_tuple>::value;

  struct ode_system {
    equation_type equation{};

    void operator()(span<double, Size + 1> dTdt_and_dcdt,
                    span<const double, Size + 1> T_and_c, double /* t */) const
        noexcept {
      // Compute production rates by calling the reaction mechanism.
      span<double, Size> dcdt = drop<1>(dTdt_and_dcdt);
      span<const double, Size> c = drop<1>(T_and_c);
      const double T = T_and_c[0];
      const double c_sum = std::accumulate(c.begin(), c.end(), 0.0);
      const double R = equation.get_universal_gas_constant();
      const double P = c_sum * R * T;
      // Fill drhoXdt here by calling the get_production_rates function of the
      // underlying mechanism.
      equation.get_production_rates(dcdt, c, T, P);

      // Get the mole fractions by computing the mean molar masses
      span<const double, Size> M = equation.get_molar_masses();
      const double mean_M = fub::transform_reduce(M, c, double(0)) / c_sum;
      std::array<double, Size> rhoX;
      std::transform(c.begin(), c.end(), rhoX.begin(),
                     [=](double c_) { return c_ * mean_M; });
      const double rho = std::accumulate(rhoX.begin(), rhoX.end(), double(0));

      // Compute dTdt such that the internal energy stays constant!
      std::array<double, Size> h;
      equation.get_specific_enthalpies_of_formation(h, T);
      std::array<double, Size> cp;
      equation.get_specific_heat_capacities_at_constant_pressure(cp, T);
      // We set dTdt as shown in Phillips thesis here
      double dTdt = 0;
      double mean_cp = 0;
      for (int i = 0; i < Size - 1; ++i) {
        dTdt += (h[i] - R / M[i] * T) * (M[i] * dcdt[i]);
        mean_cp += cp[i] * c[i] * M[i] / rho;
      }
      dTdt /= (rho * (R / mean_M - mean_cp));
      dTdt_and_dcdt[0] = dTdt;
    }
  };

  template <typename State>
  auto retrieve_T_and_c(const equation_type& equation, const State& state) const
      noexcept {
    using namespace variables;
    std::array<double, Size + 1> T_and_c;
    T_and_c[0] = equation.get(temperature, state);
    std::array<double, Size> Y = equation.get_mass_fractions(state);
    double rho = equation.get(density, state);
    span<const double, Size> molar_masses = equation.get_molar_masses();
    for (index i = 0; i < Size; ++i) {
      T_and_c[i + 1] = rho * Y[i] / molar_masses[i];
    }
    return T_and_c;
  }

  complete_state advance_state(const equation_type& equation,
                               const complete_state& state, duration dt) const {
    using namespace variables;

    // Solve the above ODE here by calling an ode solver.
    ode_system system{equation};
    std::array<double, Size + 1> T_and_c = retrieve_T_and_c(equation, state);
    m_ode_solver.integrate(system, make_span(T_and_c), dt);

    /// Make a new state from the T and Y values.
    const double T = T_and_c[0];
    const double R = equation.get_universal_gas_constant();
    span<double, Size> c = drop<1>(make_span(T_and_c));
    const double c_sum = std::accumulate(c.begin(), c.end(), 0.0);
    const double P = c_sum * R * T;

    // Get the mole fractions by computing the mean molar masses
    span<const double, Size> M = equation.get_molar_masses();
    const double mean_M = fub::transform_reduce(M, c, double(0)) / c_sum;
    std::array<double, Size> rhoX;
    std::transform(c.begin(), c.end(), rhoX.begin(),
                   [=](double c_) { return c_ * mean_M; });
    return equation.set_momentum(equation.set_TPX(T, P, rhoX),
                                 equation.get_momentum(state));
  }

  Grid step(const Grid& grid, duration dt) const {
    using node_type = typename traits::node_type;
    Grid next{grid.equation(), grid.patch_extents()};
    kinetic_source_term_detail::advance_patch_action<Grid, kinetic_source_term>
        advance_patch{};
    for (const partition_type& partition : grid) {
      const auto& octant = traits::octant(partition);
      auto where = traits::get_location(partition);
      node_type node = traits::node(partition);
      equation_type equation = grid.equation();
      next.insert(next.end(), octant,
                  traits::dataflow_action(advance_patch, std::move(where),
                                          std::move(node), *this, equation,
                                          dt));
    }
    return next;
  }

  template <typename Coordinates, typename BoundaryCondition>
  Grid step(const Grid& grid, duration dt, const Coordinates&,
            const BoundaryCondition&) const {
    return step(grid, dt);
  }

  template <typename Coordinates, typename BoundaryCondition>
  auto get_time_step_size(const Grid& grid, const Coordinates&,
                          const BoundaryCondition&) const {
    const std::chrono::duration<double> initial{
        std::numeric_limits<double>::infinity()};
    return traits::reduce(
        grid, initial,
        [](duration x, future<duration> y) { return std::min(x, y.get()); },
        [&](const partition_type& partition) {
          kinetic_source_term_detail::get_dt_action<Grid, kinetic_source_term>
              get_dt{};
          node_type node = traits::node(partition);
          auto where = traits::get_location(partition);
          return traits::dataflow_action(get_dt, std::move(where),
                                         std::move(node), *this,
                                         grid.equation());
        });
  }

  template <typename Archive> void serialize(Archive& archive, unsigned) {
    archive& m_ode_solver;
  }

  OdeSolver m_ode_solver{};
};

namespace kinetic_source_term_detail {
template <typename Grid, typename Solver> struct advance_patch {
  using traits = grid_traits<Grid>;
  using node_type = typename traits::node_type;
  using equation_type = typename traits::equation_type;
  using patch_type = typename traits::patch_type;
  using const_patch_view = typename traits::const_patch_view_type;
  using location_type = typename traits::location_type;
  using duration = std::chrono::duration<double>;
  template <typename T> using future = typename traits::template future<T>;

  static node_type invoke(node_type node, Solver solver, equation_type equation,
                          duration dt) {
    future<const_patch_view> view = node.get_patch_view();
    node_type result = traits::dataflow(
        [node, solver, equation, dt](future<const_patch_view> future_view) {
          const_patch_view view = future_view.get();
          patch_type result(view.extents());
          patch_view_t<patch_type&> dest = make_view(result);
          std::transform(view.begin(), view.end(), dest.begin(),
                         [&](const auto& state) {
                           return solver.advance_state(equation, state, dt);
                         });
          return node_type(node.get_location(), std::move(result));
        },
        std::move(view));
    return result;
  }
};

template <typename Grid, typename Solver>
struct advance_patch_action
    : grid_traits<Grid>::template make_action<
          decltype(&advance_patch<Grid, Solver>::invoke),
          &advance_patch<Grid, Solver>::invoke,
          advance_patch_action<Grid, Solver>> {};

template <typename Grid, typename Solver> struct get_dt {
  using traits = grid_traits<Grid>;
  using node_type = typename traits::node_type;
  using equation_type = typename traits::equation_type;
  using patch_type = typename traits::patch_type;
  using const_patch_view = typename traits::const_patch_view_type;
  using location_type = typename traits::location_type;
  using duration = std::chrono::duration<double>;
  template <typename T> using future = typename traits::template future<T>;

  static duration invoke(node_type node, Solver solver,
                         equation_type equation) {
    const_patch_view view = node.get_patch_view().get();
    auto get_max_dt = [=](auto& state) {
      typename kinetic_source_term<Grid>::ode_system system{equation};
      auto T_and_c = solver.retrieve_T_and_c(equation, state);
      decltype(T_and_c) dTdt_and_dcdt;
      /// Compute derivatives and retrieve them
      system(dTdt_and_dcdt, T_and_c, 0);
      auto c = drop<1>(make_span(T_and_c));
      auto dcdt = drop<1>(make_span(dTdt_and_dcdt));
      const double T = T_and_c[0];
      const double dTdt = dTdt_and_dcdt[0];
      /// Compute max dt for each component and take the minimum
      // auto molar_masses = solver.equation.get_molar_masses();
      // const double rho = solver.equation.get(Density(), state);
      const double max_dt = fub::transform_reduce(
          dcdt.begin(), dcdt.end(), c.begin(),
          std::numeric_limits<double>::infinity(),
          [](double t1, double t2) { return std::min(t1, t2); },
          [=](double dc_dt, double c) {
            return dc_dt < 0 ? std::abs(c / dc_dt)
                             : std::numeric_limits<double>::infinity();
          });
      return dTdt < 0 ? std::min(max_dt, std::abs(T / dTdt)) : max_dt;
    };
    double max_dt = std::accumulate(view.begin(), view.end(),
                                    std::numeric_limits<double>::infinity(),
                                    [=](double dt, const auto& state) {
                                      return std::min(dt, get_max_dt(state));
                                    });
    assert(max_dt > 0);
    return std::chrono::duration<double>(max_dt);
  }
};

template <typename Grid, typename Solver>
struct get_dt_action
    : grid_traits<Grid>::template make_action<
          decltype(&get_dt<Grid, Solver>::invoke),
          &get_dt<Grid, Solver>::invoke, get_dt_action<Grid, Solver>> {};
} // namespace kinetic_source_term_detail

} // namespace euler
} // namespace fub

#endif // !KINETICSOURCETERM_HPP
