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

#ifndef FUB_SOLVER_STRANG_SPLITTING_HPP
#define FUB_SOLVER_STRANG_SPLITTING_HPP

#include <chrono>

namespace fub {

struct strang_splitting {
  /// @brief This the the root case where only one solver is specified. After
  /// calling that one we return the function.
  template <typename Grid, typename Coordinates, typename BoundaryCondition,
            typename LastSolver>
  static auto step(Grid&& grid, std::chrono::duration<double> dt,
                   const Coordinates& coordinates,
                   const BoundaryCondition& boundary_condition,
                   const LastSolver& last) {
    return last.step(std::forward<Grid>(grid), dt, coordinates,
                     boundary_condition);
  }

  /// @brief This function recursively invokes each specified solver.
  template <typename Grid, typename Coordinates, typename BoundaryCondition,
            typename FirstSolver, typename... Solvers>
  static auto step(Grid&& grid, std::chrono::duration<double> dt,
                   const Coordinates& coordinates,
                   const BoundaryCondition& boundary_condition,
                   const FirstSolver& first, const Solvers&... solver) {
    auto next = first.step(std::forward<Grid>(grid), dt / 2, coordinates,
                           boundary_condition);
    next =
        step(std::move(next), dt, coordinates, boundary_condition, solver...);
    return first.step(std::move(next), dt / 2, coordinates, boundary_condition);
  }
};

} // namespace fub

#endif // !STRANG_SPLITTING_HPP
