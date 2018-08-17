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

#include "fub/equation.hpp"
#include "fub/face.hpp"
#include "fub/grid.hpp"
#include "fub/permutate_dimensions.hpp"
#include "fub/variable_view.hpp"

#include <chrono>
#include <numeric>

namespace fub {
inline namespace v1 {

template <typename Grid> struct neighbors_side_data_t {
  quadrant_t<Grid> quadrant;
  patch_buffer_t<Grid> data;
};

template <typename Grid> struct neighbors_data_t {
  neighbors_side_data_t<Grid> left;
  neighbors_side_data_t<Grid> right;
};

/// The hyperbolic system solver performs a dimensional split time step scheme.
template <typename FluxMethod, typename TimeIntegrator, typename Equation,
          typename BoundaryCondition>
struct simple_hyperbolic_system_solver {
  FluxMethod flux_method;
  TimeIntegrator time_integrator;
  Equation equation;
  BoundaryCondition boundary_condition;

  using floating_point = floating_point_t<Equation>;
  using duration = std::chrono::duration<floating_point>;

  /// Returns neighbor data after applying a boundary condition or a
  /// linear interpolation scheme.
  template <typename Grid, typename Quadrant>
  neighbors_data_t<Grid> compute_neighbor_data(const Grid& grid, Quadrant q,
                                               int dim);

  /// Estimates a time step size which does not violate the CFL condition.
  template <typename Grid>
  duration estimate_stable_time_step_size(const Grid& grid) const;

  /// Performs a conservative time step with specified time step size.
  template <typename Grid>
  Grid time_step(const Grid& grid, duration time_step_size) const;
};

template <typename FM, typename TI, typename Equation, typename BC>
simple_hyperbolic_system_solver<FM, TI, Equation, BC>
make_hyperbolic_system_solver(FM flux_method, TI time_integrator,
                              Equation equation, BC boundary_condition) {
  return {flux_method, time_integrator, equation, boundary_condition};
}

////////////////////////////////////////////////////////////////////////////////
//                                                              Implementation

template <typename Equation, typename Neighbors>
void _permutate_neighbors(Equation eq, Neighbors& neighbors,
                          std::array<int, Equation::rank()> permutation,
                          patch_t<Equation> scratch) noexcept {
  eq.permutate_dimensions(neighbors.left.data, permutation, scratch);
  std::memcpy(neighbors.left.data.span().data(), scratch.span().data(),
              scratch.span().byte_size());
  eq.permutate_dimensions(neighbors.right.data, permutation, scratch);
  std::memcpy(neighbors.right.data.span().data(), scratch.span().data(),
              scratch.span().byte_size());
}

template <typename FM, typename TI, typename Eq, typename BC>
template <typename Grid, typename Quadrant>
neighbors_data_t<Grid>
simple_hyperbolic_system_solver<FM, TI, Eq, BC>::compute_neighbor_data(
    const Grid& grid, Quadrant q, int dim) {
  neighbors_data_t<Grid> neighbors;
  // Compute Left Neighbor
  neighbors.left.data =
      patch_buffer_t<Grid>(grid.variable_list(), grid.patch_extents());
  fub::face fL{axis(dim), direction::left};
  neighbors.left.quadrant = face_neighbor(q, fL);
  auto lefts = grid.face_neighbors(q, fL);
  if (lefts.size() == 0) {
    boundary_condition.compute_neighbor_data(equation, grid, q, fL,
                                             neighbors.left.data);
  } else {
    interpolate_neighbor_data(grid, q, fL, neighbors.left.data.span());
  }
  // Compute Right Neighbor
  neighbors.right.data =
      patch_buffer_t<Grid>(grid.variable_list(), grid.patch_extents());
  fub::face fR{axis(dim), direction::right};
  neighbors.right.quadrant = face_neighbor(q, fR);
  auto rights = grid.face_neighbors(q, fR);
  if (rights.size() == 0) {
    boundary_condition.compute_neighbor_data(equation, grid, q, fR,
                                             neighbors.right.data);
  } else {
    interpolate_neighbor_data(grid, q, fR, neighbors.right.data.span());
  }
  return neighbors;
}

template <typename Solver, typename Grid, typename Quadrant>
typename Solver::duration
_estimate_stable_time_step_size_split(Solver solver, const Grid& grid,
                                      Quadrant quadrant, int dim) {
  patch_data_t<const Grid> middle = grid.patch_data(quadrant);
  patch_buffer_t<Grid> scratch(middle.get_variable_list(),
                               middle.get_extents());
  const coordinates_t<Grid> coordinates = grid.coordinates(quadrant);
  neighbors_data_t<Grid> neighbors =
      solver.compute_neighbor_data(grid, quadrant, dim);
  _permutate_neighbors(solver.equation, neighbors, {0, dim}, scratch);
  solver.equation.permutate_dimensions(middle, {0, dim}, scratch);
  middle = patch_data_t<const Grid>(scratch);
  const floating_point_t<Grid> dx = coordinates.dx(dim);
  return solver.flux_method.estimate_stable_time_step_size(
      solver.equation, dx, neighbors.left.data, middle, neighbors.right.data);
}

template <typename Solver, typename Grid, typename Quadrant>
typename Solver::duration _estimate_stable_time_step_size(Solver solver,
                                                          const Grid& grid,
                                                          Quadrant quadrant) {
  typename Solver::duration dt =
      _estimate_stable_time_step_size_split(solver, grid, quadrant, 0);
  for (int dim = 1; dim < Grid::rank(); ++dim) {
    dt = std::min(
        dt, _estimate_stable_time_step_size_split(solver, grid, quadrant, dim));
  }
  return dt;
}

template <typename FM, typename TI, typename Eq, typename BC>
template <typename Grid>
std::chrono::duration<floating_point_t<Eq>>
simple_hyperbolic_system_solver<FM, TI, Eq, BC>::estimate_stable_time_step_size(
    const Grid& grid) const {
  auto&& trees = grid.get_forest().trees();
  constexpr floating_point_t<Eq> inf =
      std::numeric_limits<floating_point_t<Eq>>::infinity();
  using duration_t = std::chrono::duration<floating_point_t<Eq>>;
  grid.exchange_ghost_data();
  const duration_t dt = std::accumulate(
      trees.begin(), trees.end(), duration_t(inf),
      [&](duration_t dt, auto&& tree) {
        auto&& quadrants = tree.quadrants();
        return std::accumulate(
            quadrants.begin(), quadrants.end(), dt,
            [&](duration_t dt, auto quad) {
              return std::min(
                  dt, _estimate_stable_time_step_size(*this, grid, quad));
            });
      });
  return dt;
}

template <typename FM, typename TI, typename Eq, typename BC, typename Grid>
Grid _time_step_split(simple_hyperbolic_system_solver<FM, TI, Eq, BC> solver,
                      const Grid& current,
                      std::chrono::duration<floating_point_t<Eq>> dt, int dim) {
  Grid next = current.clone();
  current.exchange_ghost_data();
  for_each_patch_box(next, [&](auto&& box) {
    const floating_point_t<Eq> dx = current.coordinates(box).dx(dim);
    patch_data_t<Grid> result = next.patch_data(box);
    patch_data_t<const Grid> mid = current.patch_data(box);
    neighbors_data_t<Grid> neighbors =
        solver.compute_neighbor_data(current, box, dim);
    patch_data_t<const Grid> left = neighbors.left.data;
    patch_data_t<const Grid> right = neighbors.right.data;
    patch_buffer_t<Grid> scratch;
    if (dim != 0) {
      scratch =
          patch_buffer_t<Grid>(mid.get_variable_list(), mid.get_extents());
      _permutate_neighbors(solver.equation, neighbors, {0, dim}, scratch);
      solver.equation.permutate_dimensions(mid, {0, dim}, scratch);
      mid = patch_data_t<const Grid>(scratch);
    }
    constexpr int Rank = Eq::rank();
    patch_buffer_t<Grid, conservative_t<Eq>> fluxes(
        grow(mid.get_extents(), int_c<0>));
    solver.time_integrator.integrate(solver.equation, solver.flux_method, dt,
                                     dx, left, mid, right, result, fluxes);
    if (dim != 0) {
      solver.equation.permutate_dimensions(result, {0, dim}, scratch);
      std::memcpy(result.span().data(), scratch.span().data(),
                  scratch.span().byte_size());
    }
  });
  return next;
}

template <typename FM, typename TI, typename Eq, typename BC>
template <typename Grid>
Grid simple_hyperbolic_system_solver<FM, TI, Eq, BC>::time_step(
    const Grid& grid, duration dt) const {
  constexpr int Rank = Eq::rank();
  const duration half_dt = dt / 2;
  Grid next = _time_step_split(*this, grid, half_dt, 0);
  for (int dim = 1; dim < Rank - 1; ++dim) {
    next = _time_step_split(*this, next, half_dt, dim);
  }
  next = _time_step_split(*this, next, dt, Rank - 1);
  for (int dim = Rank - 2; dim >= 0; --dim) {
    next = _time_step_split(*this, next, half_dt, dim);
  }
  return next;
}

} // namespace v1
} // namespace fub

#endif // !FUB_HYPERBOLIC_SYSTEM_SOLVER_HPP
