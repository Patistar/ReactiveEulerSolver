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

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#include <hpx/hpx.hpp>
#pragma clang diagnostic pop

#include "fub/grid.hpp"
#include "fub/octree.hpp"
#include "fub/patch.hpp"

namespace fub {
namespace hpx {

using ::hpx::async;
using ::hpx::future;
using ::hpx::make_ready_future;
using ::hpx::shared_future;

template <typename Eq, typename Extents> struct grid_node {
  using patch_type = patch<as_tuple_t<complete_state_t<Eq>>, Extents>;
  grid_node() = default;
  explicit grid_node(patch_type&& b) noexcept : patch{std::move(b)} {}
  explicit grid_node(const patch_type& b) noexcept : patch{b} {}
  patch_type patch;
};

template <typename Eq, typename PatchExtents> class grid {
public:
  static constexpr int rank = PatchExtents::rank;

  using equation_type = Eq;
  using extents_type = PatchExtents;
  using patch_type = patch<as_tuple_t<complete_state_t<Eq>>, extents_type>;
  using node_type = grid_node<equation_type, extents_type>;
  using mapped_type = shared_future<node_type>;
  using octree_type = octree<mapped_type, rank>;
  using octant_type = octant<rank>;
  using partition_type = typename octree_type::value_type;
  using iterator = typename octree_type::iterator;
  using const_iterator = typename octree_type::const_iterator;

private:
  octree_type m_tree;
  extents_type m_patch_extents;

public:
  grid() = default;

  grid(int depth, const extents_type& extents = extents_type())
      : m_tree{}, m_patch_extents{extents} {
    octant_type octant{depth, {0}};
    octant_type last = octant_type{}.upper_descendant(depth);
    while (octant != last) {
      auto data = hpx::make_ready_future(node_type{});
      m_tree.insert(m_tree.end(), {octant, std::move(data)});
      octant = octant.next();
    }
    auto pointer = hpx::make_ready_future(node_type{});
    m_tree.insert(m_tree.end(), {octant, std::move(pointer)});
  }

  // ITERATORS

  const_iterator begin() const noexcept { return m_tree.begin(); }
  const_iterator end() const noexcept { return m_tree.end(); }

  iterator begin() noexcept { return m_tree.begin(); }
  iterator end() noexcept { return m_tree.end(); }

  /// @brief Returns the partition located at the specified octant.
  /// present it returns an empty optional.
  template <axis Axis, direction Direction>
  optional<partition_type> get_face_neighbor(const octant_type& octant) const
      noexcept {
    optional<octant_type> neighbor = face_neighbor(octant, {Axis, Direction});
    if (neighbor) {
      const_iterator it = m_tree.find(*neighbor);
      if (it != m_tree.end()) {
        return *it;
      }
    }
    return nullopt;
  }

  /// @brief Returns the amount of patches which are contained in this grid.
  auto size() const noexcept { return m_tree.size(); }

  /// @brief Returns the extents of a patch contained by this grid.
  const extents_type& patch_extents() const noexcept { return m_patch_extents; }

  auto insert(const_iterator hint, const octant_type& octant,
              future<patch_type> data) {
    shared_future<node_type> node = data.then([](future<patch_type> patch) {
      return node_type(std::move(patch.get()));
    });
    return m_tree.insert(hint, std::make_pair(octant, std::move(node)));
  }

  auto insert(const octant_type& octant, future<patch_type> data) {
    shared_future<node_type> node = data.then([](future<patch_type> patch) {
      return node_type(std::move(patch.get()));
    });
    return m_tree.insert(std::make_pair(octant, std::move(node)));
  }

  template <typename... Args>
  auto emplace(const octant_type& octant, Args&&... args) {
    auto make_node = [](auto&&... args) {
      return node_type(patch_type(std::forward<decltype(args)>(args)...));
    };
    shared_future<node_type> node =
        ::hpx::async(make_node, std::forward<Args>(args)...);
    return m_tree.insert(octant, std::move(node));
  }
};

template <typename T> struct filter_patch_fn {
  decltype(auto) operator()(T&& obj) const noexcept {
    return std::forward<T>(obj);
  }
};

template <typename Eq, typename Ex> struct filter_patch_fn<grid_node<Eq, Ex>> {
  const auto& operator()(const grid_node<Eq, Ex>& node) const noexcept {
    return node.patch;
  }
  auto& operator()(grid_node<Eq, Ex>& node) const noexcept {
    return node.patch;
  }
};

template <typename Eq, typename Ex>
struct filter_patch_fn<shared_future<grid_node<Eq, Ex>>> {
  auto operator()(const shared_future<grid_node<Eq, Ex>>& mapped) const {
    return mapped.then([](const auto& m) { return filter_patch(m.get()); });
  }
};

template <int Rank, typename Eq, typename Ex>
struct filter_patch_fn<
    std::pair<const octant<Rank>, shared_future<grid_node<Eq, Ex>>>> {
  auto operator()(const typename grid<Eq, Ex>::partition_type& p) const {
    return filter_patch_fn<typename grid<Eq, Ex>::mapped_type>{}(p.second);
  }
};

template <typename T> decltype(auto) filter_patch(T&& obj) {
  return filter_patch_fn<std::decay_t<T>>{}(std::forward<T>(obj));
}

} // namespace hpx

template <typename Eq, typename Ex> struct grid_traits<hpx::grid<Eq, Ex>> {
  using grid_type = hpx::grid<Eq, Ex>;
  using extents_type = typename grid_type::extents_type;
  using partition_type = typename grid_type::partition_type;
  using octant_type = typename grid_type::octant_type;
  using mapped_type = typename grid_type::mapped_type;
  using patch_type = typename grid_type::patch_type;
  using patch_view_type = decltype(make_view(std::declval<patch_type&>()));
  using node_type = typename grid_type::node_type;

  static constexpr int rank = grid_type::rank;

  /// @brief This function attempts to transform its argument to a patch or just
  /// forwards it otherwise.
  ///
  /// @detail The following transformations will be done
  ///
  ///              Partition -> future<patch_type>
  ///                 Mapped -> future<patch_type>
  ///                   Node -> patch_type
  /// otherwise:           T -> T
  template <typename F, typename... Args>
  static decltype(auto) dataflow(F f, Args&&... args) {
    return ::hpx::dataflow(::hpx::launch::async, ::hpx::util::unwrapping(f),
                           hpx::filter_patch(std::forward<Args>(args))...);
  }

  /// @brief Returns the octant from a specified value
  static const octant_type& octant(const partition_type& partition) noexcept {
    return partition.first;
  }

  /// @brief Returns the octant for a specified parition
  static const octant_type& octant(partition_type& partition) noexcept {
    return partition.first;
  }

  /// @brief This function specialises the reduce algorithm for the grid at
  /// hand.
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
    return ::hpx::when_all(projected).then(
        [=, binary_op = std::move(binary_op)](auto ps) {
          auto&& projected = ps.get();
          return std::accumulate(
              projected.begin(), projected.end(), initial,
              [binary_op = std::move(binary_op)](auto&& x, auto&& y) {
                return fub::invoke(binary_op, x, y.get());
              });
        });
  }
};

} // namespace fub

#endif // !GRID_HPP
