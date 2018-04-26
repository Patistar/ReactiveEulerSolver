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

#ifndef FUB_HYPERBOLICSYSTEMSOURCESOLVER_HPP
#define FUB_HYPERBOLICSYSTEMSOURCESOLVER_HPP

#include <chrono>
#include <utility>

namespace fub {

template <typename Splitting, typename Hyperbolic, typename Source>
struct HyperbolicSystemSourceSolver {
  Splitting m_splitting;
  Hyperbolic m_hyperbolic_solver;
  Source m_source_solver;

  template <typename Grid, typename Coordinates, typename BoundaryCondition>
  auto get_next_time_step(const Grid& grid, double cfl,
                          std::chrono::duration<double> limited_dt,
                          const Coordinates& coordinates,
                          const BoundaryCondition& boundary_condition) const {
    auto stable_dt_adv = m_hyperbolic_solver.get_time_step_size(
        grid, coordinates, boundary_condition);
    auto stable_dt_kin = m_source_solver.get_time_step_size(grid, coordinates,
                                                            boundary_condition);
    auto advance_grid =
        [=, solver = *this](std::chrono::duration<double> stable_dt_adv,
                            std::chrono::duration<double> stable_dt_kin) {
          std::chrono::duration<double> actual_dt =
              std::min({limited_dt, cfl * stable_dt_adv, stable_dt_kin});
          auto next = solver.m_splitting.step(
              grid, actual_dt, coordinates, boundary_condition,
              solver.m_source_solver, solver.m_hyperbolic_solver);
          return std::make_pair(std::move(next), actual_dt);
        };
    return grid_traits<Grid>::dataflow(advance_grid, std::move(stable_dt_adv),
                                       std::move(stable_dt_kin));
  }
};

template <typename Split, typename Adv, typename Source>
constexpr auto make_hyperbolic_system_source_solver(Split split, Adv adv,
                                                    Source source) {
  return HyperbolicSystemSourceSolver<Split, Adv, Source>{
      std::move(split), std::move(adv), std::move(source)};
}

} // namespace fub

#endif // !HYPERBOLICSYSTEMSOURCESOLVER_HPP
