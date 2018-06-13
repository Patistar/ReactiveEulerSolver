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

#include "fub/face.hpp"
#include "fub/grid.hpp"
#include "fub/optional.hpp"
#include "fub/patch.hpp"
#include "fub/patch_view.hpp"
#include "fub/variables.hpp"
#include "fub/variant.hpp"

#include <chrono>

namespace fub {

template <typename Equation, typename FluxMethod, typename TimeIntegrator>
struct hyperbolic_system_solver {
  Equation equation;
  FluxMethod flux_method;
  TimeIntegrator integrator;

  //////////////////////////////////////////////////////////////////////////////
  // Helper Type Functions

  template <axis Dim, direction Dir, typename BC, typename G, typename C>
  using get_face_neighbor_t = decltype(
      std::declval<const BC &>().template get_face_neighbor<2, Dim, Dir>(
          std::declval<const typename grid_traits<G>::partition_type &>(),
          std::declval<const G &>(), std::declval<const C &>()));

  //////////////////////////////////////////////////////////////////////////////
  // Get Neighbor Parittion

  template <axis Axis, direction Dir, typename Grid, typename BoundaryCondition,
            typename Coordinates>
  variant<typename grid_traits<Grid>::partition_type,
          get_face_neighbor_t<Axis, Dir, BoundaryCondition, Grid, Coordinates>>
  get_face_neighbor(const typename grid_traits<Grid>::partition_type &partition,
                    const Grid &grid,
                    const BoundaryCondition &boundary_condition,
                    const Coordinates &coordinates) const {
    using variant_t = variant<
        typename grid_traits<Grid>::partition_type,
        get_face_neighbor_t<Axis, Dir, BoundaryCondition, Grid, Coordinates>>;
    const auto &octant = grid_traits<Grid>::octant(partition);
    optional<typename grid_traits<Grid>::partition_type> existing_part(
        grid.template get_face_neighbor<Axis, Dir>(octant));
    if (existing_part) {
      return variant_t{in_place_index<0>, *existing_part};
    }
    return variant_t{
        in_place_index<1>,
        boundary_condition.template get_face_neighbor<2, Axis, Dir>(
            partition, grid, coordinates)};
  }

  //////////////////////////////////////////////////////////////////////////////
  // Integrate in Time

  template <axis Axis, typename Grid, typename Coordinates,
            typename BoundaryCondition>
  Grid step_split(const Grid &current, std::chrono::duration<double> dt,
                  const Coordinates &coordinates,
                  const BoundaryCondition &boundary_condition) const {
    Grid next{};
    for (const typename grid_traits<Grid>::partition_type &partition :
         current) {
      // Find neighbor partitions for ghost cells.
      auto left_{get_face_neighbor<Axis, direction::left>(
          partition, current, boundary_condition, coordinates)};
      auto right_{get_face_neighbor<Axis, direction::right>(
          partition, current, boundary_condition, coordinates)};
      // Get the octant of the parition which gets updates.
      const auto &octant = grid_traits<Grid>::octant(partition);

      // Whatever type left and right have, we do the same thing.
      // visit dispatches `{variant<T, S>, variant<T, S>}` into for cases
      // ({T, T}, {T, S}, {S, T}, {S, S}).
      fub::visit(
          [&](auto left, auto right) mutable {
            // Create coordinates which span the specified octant.
            const auto coords = adapt(coordinates, octant);
            // Create a function which will be possibly invoked delayed.
            auto integrate_partition = [=, solver = *this](const auto &left,
                                                           const auto &middle,
                                                           const auto &right) {
              typename grid_traits<Grid>::patch_type dest(middle.extents());
              solver.integrator.template integrate<Axis>(
                  make_view(dest), make_view(left), make_view(middle),
                  make_view(right), dt, coords, solver.flux_method,
                  solver.equation);
              return dest;
            };
            // Whenever `left`, `partition` and `right` are ready, invoke
            // integrate_partition and put the result into `next`.
            next.insert(octant,
                        grid_traits<Grid>::dataflow(
                            std::move(integrate_partition), std::move(left),
                            std::move(partition), std::move(right)));
          },
          std::move(left_), std::move(right_));
    }
    return next;
  }

  template <typename Grid, typename Coordinates, typename BoundaryCondition,
            int X, int... Dims>
  Grid step(const Grid &grid, std::chrono::duration<double> dt,
            const Coordinates &coordinates,
            const BoundaryCondition &boundary_condition,
            std::integer_sequence<int, X, Dims...>) const {
    Grid next = step_split<axis(X)>(grid, dt, coordinates, boundary_condition);
    (void)std::initializer_list<int>{
        ((void)(next = step_split<axis(Dims)>(next, dt, coordinates,
                                              boundary_condition)),
         42)...};
    return next;
  }

  template <typename Grid, typename Coordinates, typename BoundaryCondition>
  Grid step(const Grid &grid, std::chrono::duration<double> dt,
            const Coordinates &coordinates,
            const BoundaryCondition &boundary_condition) const {
    static constexpr int Rank = grid_traits<Grid>::rank;
    return step(grid, dt, coordinates, boundary_condition,
                make_int_sequence<Rank>());
  }

  //////////////////////////////////////////////////////////////////////////////
  // Get Time Step Size

  template <typename Grid, typename Coordinates, typename BoundaryCondition>
  auto get_time_step_size(const Grid &current, const Coordinates &coordinates,
                          const BoundaryCondition &boundary_condition) const {
    std::chrono::duration<double> initial{
        std::numeric_limits<double>::infinity()};
    return grid_traits<Grid>::reduce(
        current, initial, [](auto x, auto y) { return std::min(x, y); },
        [&](const typename grid_traits<Grid>::partition_type &partition) {
          auto left_{get_face_neighbor<axis::x, direction::left>(
              partition, current, boundary_condition, coordinates)};
          auto right_{get_face_neighbor<axis::x, direction::right>(
              partition, current, boundary_condition, coordinates)};
          return fub::visit(
              [&](auto &&left, auto &&right) {
                const auto &octant = grid_traits<Grid>::octant(partition);
                const auto coords = adapt(coordinates, octant);
                auto get_dt = [=, solver = *this](const auto &left,
                                                  const auto &middle,
                                                  const auto &right) {
                  std::chrono::duration<double> dt =
                      solver.flux_method.get_stable_time_step(
                          equation, make_view(left), make_view(middle),
                          make_view(right), coords);
                  assert(dt.count() > 0);
                  return dt;
                };
                return grid_traits<Grid>::dataflow(
                    std::move(get_dt), std::move(left), std::move(partition),
                    std::move(right));
              },
              std::move(left_), std::move(right_));
        });
  }

  //////////////////////////////////////////////////////////////////////////////
  // Get Next Time Step

  template <typename Grid, typename Coordinates, typename BoundaryCondition>
  auto get_next_time_step(const Grid &current, double cfl,
                          const Coordinates &coordinates,
                          const BoundaryCondition &boundary_condition) const {
    auto stable_dt =
        get_time_step_size(current, coordinates, boundary_condition);
    auto do_step = [=](std::chrono::duration<double> stable_dt) {
      std::chrono::duration<double> actual_dt = cfl * stable_dt;
      return std::make_pair(
          step(current, actual_dt, coordinates, boundary_condition), actual_dt);
    };
    return grid_traits<Grid>::dataflow(do_step, std::move(stable_dt));
  }

  // Limited Version
  template <typename Grid, typename Coordinates, typename BoundaryCondition>
  auto get_next_time_step(const Grid &current, double cfl,
                          std::chrono::duration<double> limited_dt,
                          const Coordinates &coordinates,
                          const BoundaryCondition &boundary_condition) const {
    auto stable_dt =
        get_time_step_size(current, coordinates, boundary_condition);
    auto do_step = [=](std::chrono::duration<double> stable_dt) {
      std::chrono::duration<double> actual_dt =
          std::min(limited_dt, cfl * stable_dt);
      return std::make_pair(
          step(current, actual_dt, coordinates, boundary_condition), actual_dt);
    };
    return grid_traits<Grid>::dataflow(do_step, std::move(stable_dt));
  }
};

} // namespace fub

#endif // !SYSTEMSOLVER_HP
