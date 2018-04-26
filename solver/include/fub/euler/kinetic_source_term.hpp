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
#include "fub/grid.hpp"
#include "fub/ode_solver/radau.hpp"
#include "fub/patch_view.hpp"
#include "fub/variables.hpp"

#include <array>
#include <chrono>
#include <numeric>
#include <vector>

namespace fub {
namespace euler {

template <typename Mechanism, typename OdeSolver = ode_solver::Radau>
struct kinetic_source_term {
  ideal_gas<Mechanism> equation{};
  OdeSolver m_ode_solver{};

  static constexpr int Size =
      std::tuple_size<typename ideal_gas<Mechanism>::species_tuple>::value;

  struct ode_system {
    ideal_gas<Mechanism> equation{};

    void operator()(span<double, Size + 1> dTdt_and_drhoXdt,
                    span<const double, Size + 1> T_and_rhoX,
                    double /* t */) const noexcept {
      // Compute production rates by calling the reaction mechanism.
      span<double, Size> drhoXdt = drop<1>(dTdt_and_drhoXdt);
      span<const double, Size> rhoX = drop<1>(T_and_rhoX);
      const double T = T_and_rhoX[0];
      assert(T >= 0);
      const double fractions_sum =
          std::accumulate(rhoX.begin(), rhoX.end(), 0.0);
      const double R = equation.get_universal_gas_constant();
      const double P = fractions_sum * R * T;
      // Fill drhoXdt here by calling the get_production_rates function of the
      // underlying mechanism.
      equation.get_production_rates(drhoXdt, rhoX, T, P);

      // Compute dTdt such that the internal energy stays constant!
      std::array<double, Size> h;
      std::array<double, Size> cp;
      typename ideal_gas<Mechanism>::complete_state q =
          equation.set_TPX(h, cp, T, P, rhoX);
      const double mean_cp =
          equation.get_mean_specific_heat_capacity_at_constant_pressure(q);
      span<const double, Size> M = equation.get_molar_masses();
      double dTdt = 0;
      double mean_M = 0;
      for (int i = 0; i < Size - 1; ++i) {
        dTdt += (h[i] - R / M[i] * T) * (M[i] * drhoXdt[i]);
        mean_M += rhoX[i] * M[i] / fractions_sum;
      }
      // We set dTdt as shown in Phillips thesis here
      dTdt_and_drhoXdt[0] =
          dTdt / (equation.get(Density(), q) * (R / mean_M - mean_cp));
    }
  };

  template <typename State>
  auto retrieve_T_and_rhoX(const State& state) const noexcept {
    using namespace variables;
    std::array<double, Size + 1> T_and_rhoX;
    T_and_rhoX[0] = equation.get(temperature, state);
    std::array<double, Size> Y = equation.get_mass_fractions(state);
    double rho = equation.get(density, state);
    span<const double, Size> molar_masses = equation.get_molar_masses();
    const double mean_molar_mass = fub::transform_reduce(
        Y, molar_masses, double(0.), std::plus<>{}, std::divides<>{});
    for (index i = 0; i < Size; ++i) {
      const double lambda = (mean_molar_mass / molar_masses[i]);
      T_and_rhoX[i + 1] = rho * Y[i] * lambda;
    }
    return T_and_rhoX;
  }

  template <typename State>
  std::enable_if_t<
      std::is_convertible<State,
                          typename ideal_gas<Mechanism>::complete_state>::value,
      typename ideal_gas<Mechanism>::complete_state>
  advance_state(const State& state, std::chrono::duration<double> dt) const {
    using namespace variables;

    /// Solve the above ODE here by calling an ode solver.
    ode_system system{equation};
    std::array<double, Size + 1> T_and_rhoX = retrieve_T_and_rhoX(state);
    m_ode_solver.integrate(system, make_span(T_and_rhoX), dt);

    /// Make a new state from the T and Y values.
    const double T = T_and_rhoX[0];
    const double R = equation.get_universal_gas_constant();
    span<double, Size> rhoX = drop<1>(make_span(T_and_rhoX));
    const double fractions_sum = std::accumulate(rhoX.begin(), rhoX.end(), 0.0);
    const double P = fractions_sum * R * T;
    return equation.set_momentum(equation.set_TPX(T, P, rhoX),
                                 equation.get_momentum(state));
  }

  template <typename Grid>
  Grid step(const Grid& grid, std::chrono::duration<double> dt) const {
    using Patch = typename grid_traits<Grid>::patch_type;
    Grid next{};
    auto advance_data = [=, solver = *this](Patch data) {
      auto view = make_view(data);
      std::transform(
          view.begin(), view.end(), view.begin(),
          [=](const auto& state) { return solver.advance_state(state, dt); });
      return data;
    };
    for (const typename grid_traits<Grid>::partition_type& partition : grid) {
      const auto& octant = grid_traits<Grid>::octant(partition);
      next.insert(octant, grid_traits<Grid>::dataflow(advance_data, partition));
    }
    return next;
  }

  template <typename Grid, typename Coordinates, typename BoundaryCondition>
  Grid step(const Grid& grid, std::chrono::duration<double> dt,
            const Coordinates&, const BoundaryCondition&) const {
    return step(grid, dt);
  }

  template <typename Grid, typename Coordinates, typename BoundaryCondition>
  auto get_time_step_size(const Grid& grid, const Coordinates&,
                          const BoundaryCondition&) const {
    std::chrono::duration<double> initial{
        std::numeric_limits<double>::infinity()};
    return grid_traits<Grid>::reduce(
        grid, initial, [](auto x, auto y) { return std::min(x, y); },
        [&](const typename grid_traits<Grid>::partition_type& partition) {
          auto get_dt =
              [=, solver = *this](
                  const auto& patch) -> std::chrono::duration<double> {
            const auto view = make_view(patch);
            auto get_max_dt = [=](auto& state) {
              ode_system system{equation};
              auto T_and_rhoX = retrieve_T_and_rhoX(state);
              decltype(T_and_rhoX) dTdt_and_drhoXdt;
              /// Compute derivatives and retrieve them
              system(dTdt_and_drhoXdt, T_and_rhoX, 0);
              auto rhoX = drop<1>(make_span(T_and_rhoX));
              auto drhoXdt = drop<1>(make_span(dTdt_and_drhoXdt));
              const double T = T_and_rhoX[0];
              const double dTdt = dTdt_and_drhoXdt[0];
              /// Compute max dt for each component and take the minimum
              auto molar_masses = solver.equation.get_molar_masses();
              std::array<double, Size> drhoYdt;
              std::transform(drhoXdt.begin(), drhoXdt.end(),
                             molar_masses.begin(), drhoYdt.begin(),
                             [=](double x, double m) { return x * m; });
              std::array<double, Size> Y =
                  solver.equation.get_mass_fractions(state);
              const double rho = solver.equation.get(Density(), state);
              const double max_dt = fub::transform_reduce(
                  drhoYdt.begin(), drhoYdt.end(), Y.begin(),
                  std::numeric_limits<double>::infinity(),
                  [](double t1, double t2) { return std::min(t1, t2); },
                  [=](double drhoY_dt, double y) {
                    return drhoY_dt < -1E-16
                               ? std::abs(rho * y / drhoY_dt)
                               : std::numeric_limits<double>::infinity();
                  });
              return dTdt < 0 ? std::min(max_dt, std::abs(T / dTdt)) : max_dt;
            };
            const double max_dt =
                std::accumulate(view.begin(), view.end(),
                                std::numeric_limits<double>::infinity(),
                                [=](double dt, const auto& state) {
                                  return std::min(dt, get_max_dt(state));
                                });
            assert(max_dt > 0);
            return std::chrono::duration<double>(std::max(1e-8, max_dt));
          };
          return grid_traits<Grid>::dataflow(std::move(get_dt),
                                             std::move(partition));
        });
  }
};

template <typename Mechanism, typename OdeSolver = ode_solver::Radau>
constexpr auto make_kinetic_source_term(const ideal_gas<Mechanism>& eq,
                                        OdeSolver solver = OdeSolver()) {
  return kinetic_source_term<Mechanism, OdeSolver>{eq, std::move(solver)};
}

} // namespace euler
} // namespace fub

#endif // !KINETICSOURCETERM_HPP
