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

#ifndef FUB_HYPERBOLIC_SYSTEM_SOLVER_HPP
#define FUB_HYPERBOLIC_SYSTEM_SOLVER_HPP

#include <chrono>

namespace fub {

template <typename FluxMethod, typename TimeIntegrator, typename Equation>
struct hyperbolic_system_solver {
  FluxMethod flux_method;
  TimeIntegrator time_integrator;
  Equation equation;

  using fub::p4est::basic_grid;

  template <typename T> using seconds = std::chrono::duration<T>;

  /// Estimates a time step size which does not violate the CFL condition.
  template <typename Vs, int Rank, typename T>
  seconds<T>
  estimate_stable_time_step_size(const basic_grid<Vs, Rank, T>& grid) const;

  /// Perfoms a conservative time step with specified time step size in a
  /// specific split direction `dim`.
  template <typename Vs, int Rank, typename T>
  void time_step_split(const basic_grid<Vs, Rank, T>& grid, duration<T> size,
                       axis dim) const;

  /// Perfoms a conservative time step with specified time step size.
  template <typename Vs, int Rank, typename T>
  void time_step(const basic_grid<Vs, Rank, T>& grid, duration<T> size) const;
};

////////////////////////////////////////////////////////////////////////////////
//                                                              Implementation

template <typename FluxMethod, typename TimeIntegrator, typename Equation>
template <typename Vs, int Rank, typename T>
seconds<T> hyperbolic_system_solver<FluxMethod, TimeIntegrator, Equation>::
    estimate_stable_time_step_size(
        const fub::p4est::basic_grid<Vs, Rank, T>& grid) const {
  using namespace fub::p4est;
  span<const tree<Rank>> trees = grid.forest().trees();
  constexpr T inf = std::numeric_limits<T>::infinity();
  return std::accumulate(
      trees.begin(), trees.end(), seconds(inf),
      [&](seconds<T> dt, const tree<Rank>& tree) {
        span<const quadrant<Rank>> quadrants = tree.quadrants();
        return std::accumulate(quadrants.begin(), quadrants.end(), dt,
                               [&](seconds<T> dt, quadrant<Rank> quad) {
                                 return std::min(
                                     dt, detail::estimate_stable_time_step_size(
                                             *this, grid, quadrant));
                               });
      });
}

template <typename Vs, int Rank, typename T>
basic_grid<Vs, Rank, T> time_step(const basic_grid<Vs, Rank, T>& grid,
                                  duration<T> size) const {
  basic_grid<Vs, Rank, T> next(grid);
  for (int dim = 0; dim < Rank - 1; ++dim) {
    next = time_step_split(next, size / 2, axis(dim));
  }
  time_step_split(next, size, axis(Rank - 1));
  for (int dim = Rank - 2; dim >= 0; --dim) {
    next = time_step_split(next, size / 2, axis(dim));
  }
  return next;
}

} // namespace fub

#endif // !FUB_HYPERBOLIC_SYSTEM_SOLVER_HPP
