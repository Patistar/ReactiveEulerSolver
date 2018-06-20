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

#ifndef FUB_HPX_GRID_HPP
#define FUB_HPX_GRID_HPP

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#endif
#include <hpx/hpx.hpp>
#ifdef __clang__
#pragma clang diagnostic pop
#endif

#include "fub/octree.hpp"
#include "fub/patch.hpp"
#include "fub/serial/grid.hpp"
#include "fub/type_traits.hpp"

namespace fub {
namespace parallel {
using serial::dummy_location;

template <typename Eq, typename Ex> struct grid_node_patch {
  using type = patch<as_tuple_t<complete_state_t<Eq>>, Ex>;
};

template <typename Eq, typename Extents> class grid_node {
public:
  using patch_type = typename grid_node_patch<Eq, Extents>::type;
  using patch_view = decltype(make_view(std::declval<const patch_type>()));

  grid_node() = default;

  grid_node(hpx::future<grid_node> node) // NOLINT
      : m_patch{node.then(
            [](hpx::future<grid_node> n) { return n.get().m_patch; })} {}

  grid_node& operator=(hpx::future<grid_node> node) {
    m_patch =
        node.then([](hpx::future<grid_node> n) { return n.get().m_patch; });
    return *this;
  }

  explicit grid_node(dummy_location, const Extents& extents)
      : m_patch(hpx::async([](const Extents& e) { return patch_type(e); },
                           extents)) {}

  explicit grid_node(dummy_location, patch_type&& patch)
      : m_patch(hpx::make_ready_future(std::move(patch))) {}

  explicit grid_node(dummy_location, const patch_type& patch)
      : m_patch(hpx::make_ready_future(patch)) {}

  explicit grid_node(hpx::future<patch_type> patch)
      : m_patch{std::move(patch)} {}

  // Accessor for the patch data.
  hpx::future<patch_view> get_patch_view() const {
    return m_patch.then([](const hpx::shared_future<patch_type>& patch) {
      return make_view(patch.get());
    });
  }

  hpx::future<dummy_location> get_locality() const {
    return hpx::make_ready_future(dummy_location{});
  }

private:
  hpx::shared_future<patch_type> m_patch;
};

template <typename Eq, typename PatchExtents> class grid {
public:
  static constexpr int rank = PatchExtents::rank;

  using equation_type = Eq;
  using extents_type = PatchExtents;
  using patch_type = patch<as_tuple_t<complete_state_t<Eq>>, extents_type>;
  using node_type = grid_node<equation_type, extents_type>;
  using mapped_type = node_type;
  using octree_type = octree<mapped_type, rank>;
  using octant_type = octant<rank>;
  using partition_type = typename octree_type::value_type;
  using iterator = typename octree_type::iterator;
  using const_iterator = typename octree_type::const_iterator;

  grid(const equation_type& eq, const extents_type& extents) noexcept
      : m_patch_extents{extents}, m_equation{eq} {}

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
  optional<partition_type> get_face_neighbor(const octant_type& octant) const
      noexcept {
    optional<octant_type> neighbor = face_neighbor(octant, {Axis, Dir});
    if (neighbor) {
      const_iterator it = m_tree.find(*neighbor);
      if (it != m_tree.end()) {
        return *it;
      }
    }
    return nullopt;
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
  extents_type m_patch_extents{};
  equation_type m_equation{};
};

} // namespace parallel

template <typename Equation, typename Extents>
struct grid_traits<parallel::grid<Equation, Extents>> {
  using grid_type = parallel::grid<Equation, Extents>;
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
  using bind_node = parallel::grid_node<Eq, Ex>;

  using location_type = parallel::dummy_location;

  template <typename T> using future = hpx::future<T>;

  template <typename AddressType, AddressType Address, typename T>
  struct make_action {
    template <typename... Args>
    decltype(auto) operator()(Args&&... args) const {
      return fub::invoke(Address, std::forward<Args>(args)...);
    }
  };

  static constexpr int rank = grid_type::rank;

  template <typename F, typename... Args>
  static decltype(auto) dataflow(F f, Args&&... args) {
    return hpx::dataflow(hpx::launch::async, f, std::forward<Args>(args)...);
  }
  template <typename F, typename... Args>
  static decltype(auto) dataflow_action(F f, parallel::dummy_location,
                                        Args&&... args) {
    return hpx::dataflow(hpx::launch::async, f, std::forward<Args>(args)...);
  }
  template <typename F, typename... Args>
  static decltype(auto)
  dataflow_action(F f, hpx::future<parallel::dummy_location>, Args&&... args) {
    return hpx::dataflow(hpx::launch::async, f, std::forward<Args>(args)...);
  }

  /// @brief Returns the octant for a specified parition
  /// @{
  static const octant_type& octant(const partition_type& partition) noexcept {
    return partition.first;
  }

  static const octant_type& octant(partition_type& partition) noexcept {
    return partition.first;
  }
  /// @}

  static const node_type& node(const partition_type& partition) noexcept {
    return partition.second;
  }
  static node_type& node(partition_type& partition) noexcept {
    return partition.second;
  }
  static node_type&& node(partition_type&& partition) noexcept {
    return std::move(partition.second);
  }

  static hpx::future<parallel::dummy_location>
  locality(const partition_type&) noexcept {
    return hpx::make_ready_future(parallel::dummy_location{});
  }

  template <typename T, typename BinaryOp, typename Proj>
  static hpx::future<T> reduce(const grid_type& grid, T initial,
                               BinaryOp binary_op, Proj proj) {
    using Projected = invoke_result_t<Proj, const partition_type&>;
    std::vector<Projected> projected;
    projected.reserve(grid.size());
    std::transform(grid.begin(), grid.end(), std::back_inserter(projected),
                   [proj = std::move(proj)](const partition_type& partition) {
                     return fub::invoke(proj, partition);
                   });
    return hpx::when_all(projected).then(
        [=, binary_op = std::move(binary_op)](auto ps) {
          auto projected_ = ps.get();
          return std::accumulate(projected_.begin(), projected_.end(), initial,
                                 binary_op);
          // [binary_op = std::move(binary_op)](auto&& x, auto&& y) {
          //   return fub::invoke(binary_op, std::forward<decltype(x)>(x),
          //                      std::forward<decltype(y)>(y));
          // });
        });
  }
}; // namespace fub

} // namespace fub

#endif // !GRID_HPP
