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

//#include "Reactor.h"

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
  auto retrieve_T_and_c(const State& state) const noexcept {
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

  typename ideal_gas<Mechanism>::complete_state
  advance_state(const typename ideal_gas<Mechanism>::complete_state& state,
                std::chrono::duration<double> dt) const {
    using namespace variables;

    // Solve the above ODE here by calling an ode solver.
    ode_system system{equation};
    std::array<double, Size + 1> T_and_c = retrieve_T_and_c(state);
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
    // using namespace variables;
    // static FlameMasterReactor reactor{"libburke2012.dylib"};
    // reactor.setTemperature(state[temperature]);
    // reactor.setDensity(state[density]);
    // auto mass_fractions = equation.get_mass_fractions(state);
    // reactor.setMassFractions(&mass_fractions[0]);
    // reactor.advance(dt.count());
    // std::copy_n(reactor.getMassFractions(), mass_fractions.size(),
    // mass_fractions.begin());
    // return equation.set_momentum(equation.set_TPY(reactor.getTemperature(),
    // reactor.getPressure(),
    // mass_fractions),
    // equation.get_momentum(state));
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
              auto T_and_c = retrieve_T_and_c(state);
              decltype(T_and_c) dTdt_and_dcdt;
              /// Compute derivatives and retrieve them
              system(dTdt_and_dcdt, T_and_c, 0);
              auto c = drop<1>(make_span(T_and_c));
              auto dcdt = drop<1>(make_span(dTdt_and_dcdt));
              const double T = T_and_c[0];
              const double dTdt = dTdt_and_dcdt[0];
              /// Compute max dt for each component and take the minimum
              auto molar_masses = solver.equation.get_molar_masses();
              const double rho = solver.equation.get(Density(), state);
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
            double max_dt =
                std::accumulate(view.begin(), view.end(),
                                std::numeric_limits<double>::infinity(),
                                [=](double dt, const auto& state) {
                                  return std::min(dt, get_max_dt(state));
                                });
            max_dt = std::max(max_dt, 1e-7);
            assert(max_dt > 0);
            return std::chrono::duration<double>(max_dt);
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
