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

#ifndef FUB_SERIAL_GRID_HPP
#define FUB_SERIAL_GRID_HPP

#include "fub/equation.hpp"
#include "fub/grid/config.hpp"
#include "fub/octree.hpp"
#include "fub/patch.hpp"
#include "fub/patch_view.hpp"
#include "fub/serial/ready_future.hpp"

#include <memory>

namespace fub {
namespace serial {
struct dummy_location {
  dummy_location() = default;
  dummy_location(ready_future<dummy_location>) {}
};

template <typename Eq, typename E> struct grid_node_patch {
  using type = patch<as_tuple_t<complete_state_t<Eq>>, E>;
};

template <typename Eq, typename Extents> class grid_node {
public:
  using patch_type = typename grid_node_patch<Eq, Extents>::type;
  using patch_view_type = decltype(make_view(std::declval<const patch_type>()));

  grid_node(ready_future<grid_node> node)
      : m_patch{std::move(node.get().m_patch)} {}

  grid_node& operator=(ready_future<grid_node> node) {
    m_patch = std::move(node.get().m_patch);
    return *this;
  }

  grid_node(dummy_location, const Extents& extents)
      : m_patch(std::make_shared<patch_type>(extents)) {}

  grid_node(dummy_location, patch_type&& patch)
      : m_patch(std::make_shared<patch_type>(std::move(patch))) {}

  grid_node(dummy_location, const patch_type& patch)
      : m_patch(std::make_shared<patch_type>(patch)) {}

  ready_future<patch_view_type> get_patch_view() const {
    ready_future<patch_view_type> result{make_view(*m_patch)};
    return result;
  }

  dummy_location get_location() const noexcept {
    return dummy_location{};
  }

private:
  std::shared_ptr<const patch_type> m_patch;
};

template <typename Equation, typename Extents>
bool operator==(const grid_node<Equation, Extents>& n1,
                const grid_node<Equation, Extents>& n2) noexcept {
  return n1.get_patch() == n2.get_patch();
}

template <typename Equation, typename Extents>
bool operator!=(const grid_node<Equation, Extents>& n1,
                const grid_node<Equation, Extents>& n2) noexcept {
  return !(n1 == n2);
}

template <typename Equation, typename Extents> class grid {
public:
  static constexpr int rank = Extents::rank;

  using equation_type = Equation;
  using extents_type = Extents;
  using patch_type = patch<as_tuple_t<complete_state_t<Equation>>, Extents>;
  using node_type = grid_node<Equation, Extents>;
  using octant_type = fub::octant<rank>;
  using mapped_type = node_type;
  using octree_type = fub::octree<mapped_type, rank>;
  using partition_type = typename octree_type::value_type;
  using iterator = typename octree_type::iterator;
  using const_iterator = typename octree_type::const_iterator;

  // Constructors

  /// Constructs an empty grid with specified equation and patch extents.
  explicit grid(const equation_type& eq, const extents_type& extents) noexcept
      : m_equation{eq}, m_patch_extents{extents} {}

  // Iterators

  /// Returns an iterator to the first element
  /// @{
  iterator begin() noexcept { return m_tree.begin(); }
  const_iterator begin() const noexcept { return m_tree.begin(); }
  const_iterator cbegin() const noexcept { return m_tree.cbegin(); }
  /// @}

  /// Returns an iterator one past the last elemet
  /// @{
  iterator end() noexcept { return m_tree.end(); }
  const_iterator end() const noexcept { return m_tree.end(); }
  const_iterator cend() const noexcept { return m_tree.end(); }
  /// @}

  // Access Partition

  /// Returns the partition for a specified octant if found.
  template <axis Axis, direction Dir>
  const_iterator find_face_neighbor(const octant_type& octant) const
      noexcept {
    optional<octant_type> neighbor = face_neighbor(octant, {Axis, Dir});
    if (neighbor) {
      const_iterator it = m_tree.find(*neighbor);
      return it;
    }
    return cend();
  }

  // Observers

  const extents_type& patch_extents() const noexcept { return m_patch_extents; }
  const equation_type& equation() const noexcept { return m_equation; }
  auto size() const noexcept { return m_tree.size(); }

  // Modifiers

  /// Inserts a grid node which contains patch data values.
  /// @{
  iterator insert(const_iterator hint, const octant_type& octant,
                  node_type&& node) {
    return m_tree.insert(hint, std::make_pair(octant, std::move(node)));
  }

  iterator insert(const_iterator hint, const octant_type& octant,
                  const node_type& node) {
    return m_tree.insert(hint, std::make_pair(octant, std::move(node)));
  }

  std::pair<iterator, bool> insert(const octant_type& octant,
                                   node_type&& node) {
    return m_tree.insert(std::make_pair(octant, std::move(node)));
  }

  std::pair<iterator, bool> insert(const octant_type& octant,
                                   const node_type& node) {
    return m_tree.insert(std::make_pair(octant, std::move(node)));
  }
  /// @}

private:
  octree_type m_tree{};
  equation_type m_equation;
  extents_type m_patch_extents;
};
} // namespace serial

template <typename G> struct grid_traits;

template <typename Equation, typename Extents>
struct grid_traits<serial::grid<Equation, Extents>> {
  using grid_type = serial::grid<Equation, Extents>;
  using equation_type = typename grid_type::equation_type;
  using partition_type = typename grid_type::partition_type;
  using extents_type = typename grid_type::extents_type;
  using octant_type = typename grid_type::octant_type;
  using mapped_type = typename grid_type::mapped_type;
  using patch_type = typename grid_type::patch_type;
  using patch_view_type = decltype(make_view(std::declval<patch_type&>()));
  using const_patch_view_type =
      decltype(make_view(std::declval<const patch_type&>()));
  using node_type = typename grid_type::node_type;

  template <typename Eq, typename Ex>
  using bind_node = serial::grid_node<Eq, Ex>;

  using location_type = serial::dummy_location;

  template <typename T> using future = ready_future<T>;

  template <typename AddressType, AddressType Address, typename T>
  struct make_action {
    template <typename... Args>
    decltype(auto) operator()(Args&&... args) const {
      return fub::invoke(Address, std::forward<Args>(args)...);
    }
  };

  static constexpr int rank = grid_type::rank;

  template <typename F, typename... Args>
  static std::enable_if_t<is_invocable<F, Args...>::value,
                          ready_future<invoke_result_t<F, Args...>>>
  dataflow_action(F f, const serial::dummy_location&, Args&&... args) {
    return ready_future<invoke_result_t<F, Args...>>(
        fub::invoke(std::move(f), std::forward<Args>(args)...));
  }

  template <typename F, typename... Args>
  static std::enable_if_t<is_invocable<F, Args...>::value,
                          ready_future<invoke_result_t<F, Args...>>>
  dataflow_action(F f, ready_future<serial::dummy_location>, Args&&... args) {
    return ready_future<invoke_result_t<F, Args...>>(
        fub::invoke(std::move(f), std::forward<Args>(args)...));
  }

  template <typename F, typename... Args>
  static std::enable_if_t<is_invocable<F, Args...>::value,
                          ready_future<invoke_result_t<F, Args...>>>
  dataflow(F f, Args&&... args) {
    return ready_future<invoke_result_t<F, Args...>>(
        fub::invoke(std::move(f), std::forward<Args>(args)...));
  }

  /// @brief Returns the octant from a specified value
  static const octant_type& octant(const partition_type& partition) noexcept {
    return partition.first;
  }

  static const node_type& node(const partition_type& partition) noexcept {
    return partition.second;
  }
  static node_type& node(partition_type& partition) noexcept {
    return partition.second;
  }
  static node_type&& node(partition_type&& partition) noexcept {
    return std::move(partition.second);
  }

  static auto get_location(const partition_type& part) noexcept {
    return part.second.get_location();
    ;
  }

  /// @brief This function specialises the reduce algorithm for the grid at
  /// hand.
  template <typename T, typename BinaryOp, typename Proj>
  static ready_future<T> reduce(const grid_type& grid, T initial,
                                BinaryOp binary_op, Proj proj) {
    for (const partition_type& partition : grid) {
      auto&& projected = fub::invoke(proj, partition);
      initial = fub::invoke(binary_op, initial, std::move(projected));
    }
    return ready_future<T>(std::move(initial));
  }
};

} // namespace fub

#endif // !GRID_HPP
