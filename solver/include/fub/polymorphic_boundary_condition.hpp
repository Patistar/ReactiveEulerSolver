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

#ifndef FUB_POLYMORPHIC_BOUNDARY_CONDITION_HPP
#define FUB_POLYMORPHIC_BOUNDARY_CONDITION_HPP

#include "fub/face.hpp"
#include "fub/uniform_cartesian_coordinates.hpp"

#include <memory>

namespace fub {
namespace detail {

/// This type trait
template <typename Grid, int Width, axis Axis>
struct boundary_condition_result {
  using equation_type = typename grid_traits<Grid>::equation_type;
  using extents_type = typename grid_traits<Grid>::extents_type;
  using node_type = typename grid_traits<Grid>::node_type;
  using reduced_extents_type = decltype(
      replace_extent(extents_type(), int_c<as_int(Axis)>, int_c<Width>));
  using result_type =
      typename grid_traits<Grid>::template bind_node<equation_type,
                                                     reduced_extents_type>;

  using type = result_type;
};

template <typename Grid, int Width, axis Axis>
using boundary_condition_result_t =
    typename boundary_condition_result<Grid, Width, Axis>::type;

template <typename Grid, int Width, axis...> struct boundary_condition;

template <typename Grid, int Width, axis Axis>
struct boundary_condition<Grid, Width, Axis> {
  using partition_type = typename grid_traits<Grid>::partition_type;

  virtual boundary_condition_result_t<Grid, Width, Axis> get_face_neighbor(
      const partition_type& partition, const Grid& grid,
      const uniform_cartesian_coordinates<Grid::rank>& coordinates,
      direction dir, axis_constant<Axis>) const = 0;

  virtual std::unique_ptr<boundary_condition> clone() const = 0;
  virtual bool equal_to(const boundary_condition& other) const noexcept = 0;
};

template <typename Grid, int Width>
struct boundary_condition<Grid, Width, axis::x, axis::y> {
  using partition_type = typename grid_traits<Grid>::partition_type;

  virtual boundary_condition_result_t<Grid, Width, axis::x> get_face_neighbor(
      const partition_type& partition, const Grid& grid,
      const uniform_cartesian_coordinates<Grid::rank>& coordinates,
      direction dir, axis_constant<axis::x>) const = 0;

  virtual boundary_condition_result_t<Grid, Width, axis::y> get_face_neighbor(
      const partition_type& partition, const Grid& grid,
      const uniform_cartesian_coordinates<Grid::rank>& coordinates,
      direction dir, axis_constant<axis::y>) const = 0;

  virtual std::unique_ptr<boundary_condition> clone() const = 0;
  virtual bool equal_to(const boundary_condition& other) const noexcept = 0;
};

template <typename Grid, int Width>
struct boundary_condition<Grid, Width, axis::x, axis::y, axis::z> {
  using partition_type = typename grid_traits<Grid>::partition_type;

  virtual boundary_condition_result_t<Grid, Width, axis::x> get_face_neighbor(
      const partition_type& partition, const Grid& grid,
      const uniform_cartesian_coordinates<Grid::rank>& coordinates,
      direction dir, axis_constant<axis::x>) const = 0;

  virtual boundary_condition_result_t<Grid, Width, axis::y> get_face_neighbor(
      const partition_type& partition, const Grid& grid,
      const uniform_cartesian_coordinates<Grid::rank>& coordinates,
      direction dir, axis_constant<axis::y>) const = 0;

  virtual boundary_condition_result_t<Grid, Width, axis::z> get_face_neighbor(
      const partition_type& partition, const Grid& grid,
      const uniform_cartesian_coordinates<Grid::rank>& coordinates,
      direction dir, axis_constant<axis::z>) const = 0;

  virtual std::unique_ptr<boundary_condition> clone() const = 0;
  virtual bool equal_to(const boundary_condition& other) const noexcept = 0;
};

template <typename BC, typename Grid, int Width, axis... Axis>
struct boundary_condition_facade;

template <typename BC, typename Grid, int Width, axis Axis>
struct boundary_condition_facade<BC, Grid, Width, Axis>
    : boundary_condition<Grid, Width, Axis> {
  using partition_type = typename grid_traits<Grid>::partition_type;

  BC m_boundary_condition;

  boundary_condition_facade(const BC& bc) : m_boundary_condition{bc} {}
  boundary_condition_facade(BC&& bc) : m_boundary_condition{std::move(bc)} {}

  virtual boundary_condition_result_t<Grid, Width, Axis> get_face_neighbor(
      const partition_type& partition, const Grid& grid,
      const uniform_cartesian_coordinates<Grid::rank>& coordinates,
      direction dir, axis_constant<Axis>) const override {
    if (dir == direction::left) {
      return m_boundary_condition
          .template get_face_neighbor<Width, Axis, direction::left>(
              partition, grid, coordinates);
    } else {
      return m_boundary_condition
          .template get_face_neighbor<Width, Axis, direction::right>(
              partition, grid, coordinates);
    }
  }

  std::unique_ptr<boundary_condition<Grid, Width, Axis>>
  clone() const override {
    return std::make_unique<boundary_condition_facade>(m_boundary_condition);
  }

  bool equal_to(const boundary_condition<Grid, Width, Axis>& other) const
      noexcept override {
    if (auto other_ = dynamic_cast<const boundary_condition_facade*>(&other)) {
      return other_->m_boundary_condition == m_boundary_condition;
    }
    return false;
  }
};

template <typename Grid, int Width, typename BC>
struct boundary_condition_facade<BC, Grid, Width, axis::x, axis::y>
    : boundary_condition<Grid, Width, axis::x, axis::y> {
  using partition_type = typename grid_traits<Grid>::partition_type;

  BC m_boundary_condition;

  boundary_condition_facade(const BC& bc) : m_boundary_condition{bc} {}
  boundary_condition_facade(BC&& bc) : m_boundary_condition{std::move(bc)} {}

  virtual boundary_condition_result_t<Grid, Width, axis::x> get_face_neighbor(
      const partition_type& partition, const Grid& grid,
      const uniform_cartesian_coordinates<Grid::rank>& coordinates,
      direction dir, axis_constant<axis::x>) const override {
    if (dir == direction::left) {
      return m_boundary_condition
          .template get_face_neighbor<Width, axis::x, direction::left>(
              partition, grid, coordinates);
    } else {
      return m_boundary_condition
          .template get_face_neighbor<Width, axis::x, direction::right>(
              partition, grid, coordinates);
    }
  }

  virtual boundary_condition_result_t<Grid, Width, axis::y> get_face_neighbor(
      const partition_type& partition, const Grid& grid,
      const uniform_cartesian_coordinates<Grid::rank>& coordinates,
      direction dir, axis_constant<axis::y>) const override {
    if (dir == direction::left) {
      return m_boundary_condition
          .template get_face_neighbor<Width, axis::y, direction::left>(
              partition, grid, coordinates);
    } else {
      return m_boundary_condition
          .template get_face_neighbor<Width, axis::y, direction::right>(
              partition, grid, coordinates);
    }
  }

  std::unique_ptr<boundary_condition<Grid, Width, axis::x, axis::y>>
  clone() const override {
    return std::make_unique<boundary_condition_facade>(m_boundary_condition);
  }

  bool
  equal_to(const boundary_condition<Grid, Width, axis::x, axis::y>& other) const
      noexcept override {
    if (auto other_ = dynamic_cast<const boundary_condition_facade*>(&other)) {
      return other_->m_boundary_condition == m_boundary_condition;
    }
    return false;
  }
};

template <typename Grid, int Width, typename BC>
struct boundary_condition_facade<BC, Grid, Width, axis::x, axis::y, axis::z>
    : boundary_condition<Grid, Width, axis::x, axis::y> {
  using partition_type = typename grid_traits<Grid>::partition_type;

  BC m_boundary_condition;

  boundary_condition_facade(const BC& bc) : m_boundary_condition{bc} {}
  boundary_condition_facade(BC&& bc) : m_boundary_condition{std::move(bc)} {}

  virtual boundary_condition_result_t<Grid, Width, axis::x> get_face_neighbor(
      const partition_type& partition, const Grid& grid,
      const uniform_cartesian_coordinates<Grid::rank>& coordinates,
      direction dir, axis_constant<axis::x>) const override {
    if (dir == direction::left) {
      return m_boundary_condition
          .template get_face_neighbor<Width, axis::x, direction::left>(
              partition, grid, coordinates);
    } else {
      return m_boundary_condition
          .template get_face_neighbor<Width, axis::x, direction::right>(
              partition, grid, coordinates);
    }
  }

  virtual boundary_condition_result_t<Grid, Width, axis::y> get_face_neighbor(
      const partition_type& partition, const Grid& grid,
      const uniform_cartesian_coordinates<Grid::rank>& coordinates,
      direction dir, axis_constant<axis::y>) const override {
    if (dir == direction::left) {
      return m_boundary_condition
          .template get_face_neighbor<Width, axis::y, direction::left>(
              partition, grid, coordinates);
    } else {
      return m_boundary_condition
          .template get_face_neighbor<Width, axis::y, direction::right>(
              partition, grid, coordinates);
    }
  }

  virtual boundary_condition_result_t<Grid, Width, axis::z> get_face_neighbor(
      const partition_type& partition, const Grid& grid,
      const uniform_cartesian_coordinates<Grid::rank>& coordinates,
      direction dir, axis_constant<axis::z>) const override {
    if (dir == direction::left) {
      return m_boundary_condition
          .template get_face_neighbor<Width, axis::z, direction::left>(
              partition, grid, coordinates);
    } else {
      return m_boundary_condition
          .template get_face_neighbor<Width, axis::z, direction::right>(
              partition, grid, coordinates);
    }
  }

  std::unique_ptr<boundary_condition<Grid, Width, axis::x, axis::y, axis::z>>
  clone() const override {
    return std::make_unique<boundary_condition_facade>(m_boundary_condition);
  }

  bool equal_to(const boundary_condition<Grid, Width, axis::x, axis::y,
                                         axis::z>& other) const
      noexcept override {
    if (auto other_ = dynamic_cast<const boundary_condition_facade*>(&other)) {
      return other_->m_boundary_condition == m_boundary_condition;
    }
    return false;
  }
};
} // namespace detail

template <typename Grid, int Width, axis... Axis>
class polymorphic_boundary_condition;

template <typename Grid, int Width, axis Axis>
class polymorphic_boundary_condition<Grid, Width, Axis> {
private:
  using pointer =
      std::unique_ptr<detail::boundary_condition<Grid, Width, Axis>>;
  pointer m_impl;

  pointer clone() const {
    if (m_impl) {
      return m_impl->clone();
    }
    return pointer{};
  }

public:
  using partition_type = typename grid_traits<Grid>::partition_type;

  polymorphic_boundary_condition() = delete;

  polymorphic_boundary_condition(const polymorphic_boundary_condition& other)
      : m_impl{other.clone()} {}

  polymorphic_boundary_condition(polymorphic_boundary_condition&& other)
      : m_impl{std::move(other.m_impl)} {}

  template <typename BC,
            typename = std::enable_if_t<!std::is_same<
                std::decay_t<BC>, polymorphic_boundary_condition>::value>>
  polymorphic_boundary_condition(BC&& condition)
      : m_impl{std::make_unique<detail::boundary_condition_facade<
            std::decay_t<BC>, Grid, Width, Axis>>(
            std::forward<BC>(condition))} {}

  template <int W, axis A, direction Direction>
  std::enable_if_t<(W == Width && A == Axis),
                   detail::boundary_condition_result_t<Grid, Width, Axis>>
  get_face_neighbor(
      const partition_type& partition, const Grid& grid,
      const uniform_cartesian_coordinates<Grid::rank>& coordinates) const {
    if (m_impl) {
      return m_impl->get_face_neighbor(partition, grid, coordinates, Direction,
                                       axis_c<A>);
    }
    throw std::runtime_error{"polymorphic_boundary_condition::get_face_"
                             "neighbor: Empty boundary condition."};
  }

  friend bool operator==(const polymorphic_boundary_condition& b1,
                         const polymorphic_boundary_condition& b2) noexcept {
    if (b1.m_impl && b2.m_impl) {
      return b1.m_impl->equal_to(*b2.m_impl);
    } else if (!b1.m_impl && !b2.m_impl) {
      return true;
    }
    return false;
  }

  friend bool operator!=(const polymorphic_boundary_condition& b1,
                         const polymorphic_boundary_condition& b2) noexcept {
    return !(b1 == b2);
  }
};

template <typename Grid, int Width>
class polymorphic_boundary_condition<Grid, Width, axis::x, axis::y> {
private:
  using pointer = std::unique_ptr<
      detail::boundary_condition<Grid, Width, axis::x, axis::y>>;
  pointer m_impl;

  pointer clone() const {
    if (m_impl) {
      return m_impl->clone();
    }
    return pointer{};
  }

public:
  using partition_type = typename grid_traits<Grid>::partition_type;

  polymorphic_boundary_condition() = delete;

  polymorphic_boundary_condition(const polymorphic_boundary_condition& other)
      : m_impl{other.clone()} {}

  polymorphic_boundary_condition(polymorphic_boundary_condition&& other)
      : m_impl{std::move(other.m_impl)} {}

  template <typename BC,
            typename = std::enable_if_t<!std::is_same<
                std::decay_t<BC>, polymorphic_boundary_condition>::value>>
  polymorphic_boundary_condition(BC&& condition)
      : m_impl{std::make_unique<detail::boundary_condition_facade<
            std::decay_t<BC>, Grid, Width, axis::x, axis::y>>(
            std::forward<BC>(condition))} {}

  template <int W, axis Axis, direction Direction>
  std::enable_if_t<(W == Width),
                   detail::boundary_condition_result_t<Grid, Width, Axis>>
  get_face_neighbor(
      const partition_type& partition, const Grid& grid,
      const uniform_cartesian_coordinates<Grid::rank>& coordinates) const {
    if (m_impl) {
      return m_impl->get_face_neighbor(partition, grid, coordinates, Direction,
                                       axis_c<Axis>);
    }
    throw std::runtime_error{"polymorphic_boundary_condition::get_face_"
                             "neighbor: Empty boundary condition."};
  }

  friend bool operator==(const polymorphic_boundary_condition& b1,
                         const polymorphic_boundary_condition& b2) noexcept {
    if (b1.m_impl && b2.m_impl) {
      return b1.m_impl->equal_to(*b2.m_impl);
    } else if (!b1.m_impl && !b2.m_impl) {
      return true;
    }
    return false;
  }

  friend bool operator!=(const polymorphic_boundary_condition& b1,
                         const polymorphic_boundary_condition& b2) noexcept {
    return !(b1 == b2);
  }
};

} // namespace fub

#endif
