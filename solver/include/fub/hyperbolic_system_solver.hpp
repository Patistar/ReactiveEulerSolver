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
#include "fub/optional.hpp"
#include "fub/patch.hpp"
#include "fub/patch_view.hpp"
#include "fub/serial/grid.hpp"
#include "fub/variables.hpp"
#include "fub/variant.hpp"

#include <chrono>

namespace fub {
namespace detail {
/// We use this class to declare an action which integrates a patch in time.
///
/// @example TODO
///
/// @see Actions
/// @{
template <typename Grid, typename Solver, typename Coordinates, typename Left,
          typename Right, axis Axis>
struct integrate_patch;

template <typename Grid, typename Solver, typename Coordinates, typename Left,
          typename Right, axis Axis>
struct integrate_patch_action;
/// @}

/// We use this class to declare an action which gets a stable time step size
/// for patch.
///
/// @example TODO
///
/// @see Actions
/// @{
template <typename Grid, typename Solver, typename Coordinates, typename Left,
          typename Right>
struct get_dt;

template <typename Grid, typename Solver, typename Coordinates, typename Left,
          typename Right>
struct get_dt_action;
/// @}
} // namespace detail

template <typename Grid, typename BoundaryCondition, typename Coordinates,
          typename FluxMethod, typename TimeIntegrator>
struct hyperbolic_system_solver {
  // Data related types
  using traits = grid_traits<Grid>;
  static constexpr int rank = traits::rank;
  using equation_type = typename traits::equation_type;
  using partition_type = typename traits::partition_type;
  using node_type = typename traits::node_type;

  // Execution related types
  using location_type = typename traits::location_type;
  template <typename T> using future = typename traits::template future<T>;

  template <typename Left, typename Right, axis Axis>
  using integrate_patch_action =
      detail::integrate_patch_action<Grid, hyperbolic_system_solver,
                                     Coordinates, Left, Right, Axis>;

  template <typename Left, typename Right>
  using get_dt_action = detail::get_dt_action<Grid, hyperbolic_system_solver,
                                              Coordinates, Left, Right>;

  // Utility typedef
  using duration = std::chrono::duration<double>;

  //////////////////////////////////////////////////////////////////////////////
  // Helper Type Functions

  template <axis Dim, direction Dir>
  using get_face_neighbor_t = decltype(
      std::declval<const BoundaryCondition&>()
          .template get_face_neighbor<2, Dim, Dir>(
              std::declval<partition_type&>(), std::declval<const Grid&>(),
              std::declval<const Coordinates&>()));

  //////////////////////////////////////////////////////////////////////////////
  // Get Neighbor Parittion

  template <axis Axis, direction Dir>
  variant<node_type, get_face_neighbor_t<Axis, Dir>>
  get_face_neighbor(const partition_type& partition, const Grid& grid,
                    const BoundaryCondition& boundary_condition,
                    const Coordinates& coordinates) const {
    using variant_t = variant<node_type, get_face_neighbor_t<Axis, Dir>>;
    auto&& octant = traits::octant(partition);
    auto existing_part = 
        grid.template find_face_neighbor<Axis, Dir>(octant);
    if (existing_part != grid.end()) {
      return variant_t{in_place_index<0>, traits::node(*existing_part)};
    }
    return variant_t{
        in_place_index<1>,
        boundary_condition.template get_face_neighbor<2, Axis, Dir>(
            partition, grid, coordinates)};
  }

  //////////////////////////////////////////////////////////////////////////////
  // Integrate in Time

  template <axis Axis>
  Grid step_split(const Grid& current, duration dt,
                  const Coordinates& coordinates,
                  const BoundaryCondition& boundary_condition) const {
    equation_type equation = current.equation();
    Grid next{equation, current.patch_extents()};
    for (const partition_type& partition : current) {
      // Find neighbor partitions for ghost cells.
      auto left_{get_face_neighbor<Axis, direction::left>(
          partition, current, boundary_condition, coordinates)};
      auto right_{get_face_neighbor<Axis, direction::right>(
          partition, current, boundary_condition, coordinates)};
      // Fetch some Meta Data about what and where to integrate.
      const octant<rank> octant = traits::octant(partition);
      const Coordinates adapted = adapt(coordinates, octant);
      node_type node = traits::node(partition);
      location_type where = traits::locality(partition);
      // Whatever type left and right have, we do the same thing.
      // visit dispatches `{variant<T, S>, variant<T, S>}` into for cases
      // ({T, T}, {T, S}, {S, T}, {S, S}).
      fub::visit(
          [&](auto left, auto right) mutable {
            // Create a function which will be possibly invoked delayed.
            using L = std::decay_t<decltype(left)>;
            using R = std::decay_t<decltype(right)>;
            // Whenever `left`, `partition` and `right` are ready, invoke
            // integrate_partition and put the result into `next`.
            next.insert(next.end(), octant,
                        traits::dataflow_action(
                            integrate_patch_action<L, R, Axis>(), std::move(where), std::move(left),
                            std::move(node), std::move(right), *this, equation,
                            adapted, dt));
          },
          std::move(left_), std::move(right_));
    }
    return next;
  }

  template <int X, int... Dims>
  Grid step(const Grid& grid, duration dt, const Coordinates& coordinates,
            const BoundaryCondition& boundary_condition,
            std::integer_sequence<int, X, Dims...>) const {
    Grid next = step_split<axis(X)>(grid, dt, coordinates, boundary_condition);
    (void)std::initializer_list<int>{
        ((void)(next = step_split<axis(Dims)>(next, dt, coordinates,
                                              boundary_condition)),
         42)...};
    return next;
  }

  Grid step(const Grid& grid, duration dt, const Coordinates& coordinates,
            const BoundaryCondition& boundary_condition) const {
    static constexpr int Rank = traits::rank;
    return step(grid, dt, coordinates, boundary_condition,
                make_int_sequence<Rank>());
  }

  //////////////////////////////////////////////////////////////////////////////
  // Get Time Step Size

  future<duration>
  get_time_step_size(const Grid& current, const Coordinates& coordinates,
                     const BoundaryCondition& boundary_condition) const {
    duration initial{std::numeric_limits<double>::infinity()};
    return traits::reduce(
        current, initial,
        [](duration x, future<duration> y) { return std::min(x, y.get()); },
        [&](const partition_type& partition) {
          auto left_nb = get_face_neighbor<axis::x, direction::left>(
              partition, current, boundary_condition, coordinates);
          auto right_nb = get_face_neighbor<axis::x, direction::right>(
              partition, current, boundary_condition, coordinates);
          octant<rank> octant = traits::octant(partition);
          Coordinates adapted = adapt(coordinates, octant);
          node_type node = traits::node(partition);
          location_type where = traits::locality(partition);
          return fub::visit(
              [&](auto left, auto right) {
                using L = std::decay_t<decltype(left)>;
                using R = std::decay_t<decltype(right)>;
                return traits::dataflow_action(
                    get_dt_action<L, R>(), std::move(where), std::move(left),
                    std::move(node), std::move(right), *this,
                    current.equation(), adapted);
              },
              std::move(left_nb), std::move(right_nb));
        });
  }

  //////////////////////////////////////////////////////////////////////////////
  // Get Next Time Step

  future<std::pair<Grid, duration>>
  get_next_time_step(const Grid& current, double cfl,
                     const Coordinates& coordinates,
                     const BoundaryCondition& boundary_condition) const {
    auto stable_dt =
        get_time_step_size(current, coordinates, boundary_condition);
    auto do_step = [=](future<duration> stable_dt) {
      duration actual_dt = cfl * stable_dt.get();
      return std::make_pair(
          step(current, actual_dt, coordinates, boundary_condition), actual_dt);
    };
    return traits::dataflow(do_step, std::move(stable_dt));
  }

  // Limited Version

  future<std::pair<Grid, duration>>
  get_next_time_step(const Grid& current, double cfl, duration limited_dt,
                     const Coordinates& coordinates,
                     const BoundaryCondition& boundary_condition) const {
    auto stable_dt =
        get_time_step_size(current, coordinates, boundary_condition);
    auto do_step = [=](future<duration> stable_dt) {
      duration actual_dt = std::min(limited_dt, cfl * stable_dt.get());
      return std::make_pair(
          step(current, actual_dt, coordinates, boundary_condition), actual_dt);
    };
    return traits::dataflow(do_step, std::move(stable_dt));
  }

  template <typename Archive> void serialize(Archive& archive, unsigned) {
    archive & flux_method;
    archive & integrator;
  }

  FluxMethod flux_method;
  TimeIntegrator integrator;
};

namespace detail {
//////////////////////////////////////////////////////////////////////////////
// Integrate Patch Action

template <typename Grid, typename Solver, typename Coordinates, typename Left,
          typename Right, axis Axis>
struct integrate_patch {
  using traits = grid_traits<Grid>;
  static constexpr int rank = traits::rank;
  using equation_type = typename traits::equation_type;
  using node_type = typename traits::node_type;
  using duration = std::chrono::duration<double>;

  static node_type invoke(Left left, node_type middle,
                          Right right, Solver solver,
                          equation_type equation,
                          Coordinates coords, duration dt) {
    auto lv = left.get_patch_view();
    auto mv = middle.get_patch_view();
    auto rv = right.get_patch_view();
    node_type result = traits::dataflow(
        [left, middle, right, solver, equation, coords, dt](auto lv, auto mv,
                                                            auto rv) {
          auto left_view = lv.get();
          auto middle_view = mv.get();
          auto right_view = rv.get();
          typename traits::patch_type dest(middle_view.extents());
          solver.integrator.template integrate<Axis>(
              make_view(dest), left_view, middle_view, right_view, dt, coords,
              solver.flux_method, equation);
          node_type result(middle.get_locality(), std::move(dest));
          return result;
        }, std::move(lv), std::move(mv), std::move(rv));
    return result;
  };
};

template <typename Grid, typename Solver, typename Coordinates, typename Left,
          typename Right, axis Axis>
struct integrate_patch_action
    : grid_traits<Grid>::template make_action<
          decltype(&integrate_patch<Grid, Solver, Coordinates, Left, Right,
                                    Axis>::invoke),
          &integrate_patch<Grid, Solver, Coordinates, Left, Right,
                           Axis>::invoke,
          integrate_patch_action<Grid, Solver, Coordinates, Left, Right,
                                 Axis>> {};

//////////////////////////////////////////////////////////////////////////////
// Get Time Step Size Action

template <typename Grid, typename Solver, typename Coordinates, typename Left,
          typename Right>
struct get_dt {
  using traits = grid_traits<Grid>;

  static constexpr int rank = traits::rank;
  using equation_type = typename traits::equation_type;
  using node_type = typename traits::node_type;
  using duration = std::chrono::duration<double>;
  template <typename T> using future = typename traits::template future<T>;

  static future<duration> invoke(Left left, node_type middle, Right right,
                                 Solver solver, equation_type equation,
                                 Coordinates coords) {
    auto lview = left.get_patch_view();
    auto mview = middle.get_patch_view();
    auto rview = right.get_patch_view();
    return traits::dataflow(
        [left, middle, right, solver, equation, coords](auto l, auto m,
                                                        auto r) {
          duration dt = solver.flux_method.get_stable_time_step(
              equation, l.get(), m.get(), r.get(), coords);
          assert(dt.count() > 0);
          return dt;
        },
        std::move(lview), std::move(mview), std::move(rview));
  }
};

template <typename Grid, typename Solver, typename Coordinates, typename Left,
          typename Right>
struct get_dt_action
    : grid_traits<Grid>::template make_action<
          decltype(&get_dt<Grid, Solver, Coordinates, Left, Right>::invoke),
          &get_dt<Grid, Solver, Coordinates, Left, Right>::invoke,
          get_dt_action<Grid, Solver, Coordinates, Left, Right>> {};
} // namespace detail

} // namespace fub

#endif // !SYSTEMSOLVER_HP
