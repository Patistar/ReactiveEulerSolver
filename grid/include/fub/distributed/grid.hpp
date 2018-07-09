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

#ifndef FUB_DISTRIBUTED_GRID_HPP
#define FUB_DISTRIBUTED_GRID_HPP

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#endif
#include <hpx/include/components.hpp>
#include <hpx/lcos/when_all.hpp>
#ifdef __clang__
#pragma clang diagnostic pop
#endif

#include "fub/interval_map.hpp"
#include "fub/octree.hpp"
#include "fub/optional.hpp"
#include "fub/patch.hpp"
#include "fub/serial/grid.hpp"

namespace fub {
namespace distributed {
/// @brief This meta function computes the type of the data patches.
template <typename Eq, typename E> struct grid_patch {
  using type = patch<as_tuple_t<complete_state_t<Eq>>, E>;
};

/// @brief This is a server side grid node hpx component.
///
/// This type lives in the AGAS (Active Global Address Space).
template <typename Eq, typename Extents>
class grid_node_server
    : public hpx::components::component_base<grid_node_server<Eq, Extents>> {
public:
  using patch_type = typename grid_patch<Eq, Extents>::type;

  grid_node_server() = default;

  explicit grid_node_server(const Extents& extents) : m_patch(extents) {}

  explicit grid_node_server(patch_type&& patch) : m_patch(std::move(patch)) {}

  explicit grid_node_server(const patch_type& patch) : m_patch(patch) {}

  patch_type get_patch() const { return m_patch; }

  HPX_DEFINE_COMPONENT_ACTION(grid_node_server, get_patch);

private:
  patch_type m_patch;
};

#define FUB_REGISTER_GRID_NODE_COMPONENT_DECLARATION(Eq, Ext)                  \
  using grid_node_server_get_patch_action_type =                               \
      ::fub::distributed::grid_node_server<Eq, Ext>::get_patch_action;         \
  HPX_REGISTER_ACTION_DECLARATION(grid_node_server_get_patch_action_type)

#define FUB_REGISTER_GRID_NODE_COMPONENT(Eq, Ext)                              \
  using grid_node_server_get_patch_action_ghost_type =                         \
      ::fub::distributed::grid_node_server<Eq,                                 \
                                           fub::extents<2>>::get_patch_action; \
  HPX_REGISTER_ACTION(grid_node_server_get_patch_action_ghost_type);           \
  using grid_node_component_ghost_type = ::hpx::components::component<         \
      ::fub::distributed::grid_node_server<Eq, fub::extents<2>>>;              \
  HPX_REGISTER_COMPONENT(grid_node_component_ghost_type)                       \
  using grid_node_server_get_patch_action_type =                               \
      ::fub::distributed::grid_node_server<Eq, Ext>::get_patch_action;         \
  HPX_REGISTER_ACTION(grid_node_server_get_patch_action_type);                 \
  using grid_node_component_type = ::hpx::components::component<               \
      ::fub::distributed::grid_node_server<Eq, Ext>>;                          \
  HPX_REGISTER_COMPONENT(grid_node_component_type)

template <typename Eq, typename Extents>
class grid_node
    : public hpx::components::client_base<grid_node<Eq, Extents>,
                                          grid_node_server<Eq, Extents>> {
  using base_type = hpx::components::client_base<grid_node<Eq, Extents>,
                                                 grid_node_server<Eq, Extents>>;

public:
  using patch_type = typename grid_patch<Eq, Extents>::type;
  using patch_view = decltype(make_view(std::declval<const patch_type>()));

  grid_node() = default;

  /// Create new component on locality 'where' and initialize the held data
  grid_node(hpx::id_type where, const Extents& extents)
      : base_type(hpx::new_<grid_node_server<Eq, Extents>>(std::move(where),
                                                           extents)) {}

  /// Create new component on locality 'where' and initialize the held data
  grid_node(hpx::id_type where, const patch_type& patch)
      : base_type(hpx::new_<grid_node_server<Eq, Extents>>(std::move(where),
                                                           patch)) {}

  /// Create new component on locality 'where' and initialize the held data
  grid_node(hpx::id_type where, patch_type&& patch)
      : base_type(hpx::new_<grid_node_server<Eq, Extents>>(std::move(where),
                                                           std::move(patch))) {}

  /// Unwrap a future node
  grid_node(hpx::future<grid_node> node) : base_type(std::move(node)) {}
  grid_node(hpx::future<hpx::id_type>&& id) : base_type(std::move(id)) {}
  grid_node(const hpx::future<hpx::id_type>& id) : base_type(id) {}

  hpx::id_type get_location() const { return hpx::get_colocation_id(base_type::get_id()).get(); }

  void retrieve_patch() {
    if (!m_patch.valid()) {
      m_patch =
          hpx::async(typename grid_node_server<Eq, Extents>::get_patch_action(),
                     base_type::get_id());
    }
  }

  hpx::future<patch_view> get_patch_view() {
    if (!m_patch.valid()) {
      retrieve_patch();
    }
    return m_patch.then([](const hpx::shared_future<patch_type>& p) {
      return make_view(p.get());
    });
  }

private:
  hpx::shared_future<patch_type> m_patch;
};

template <typename Eq, typename PatchExtents> class grid {
public:
  static constexpr int rank = PatchExtents::rank;

  using equation_type = Eq;
  using extents_type = PatchExtents;
  using node_type = grid_node<equation_type, extents_type>;
  using patch_type = typename node_type::patch_type;
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
    return m_tree.insert(hint, std::make_pair(octant, node));
  }

  std::pair<iterator, bool> insert(const octant_type& octant,
                                   node_type&& node) {
    return m_tree.insert(std::make_pair(octant, std::move(node)));
  }

  std::pair<iterator, bool> insert(const octant_type& octant,
                                   const node_type& node) {
    return m_tree.insert(std::make_pair(octant, node));
  }
  /// @}

private:
  octree_type m_tree{};
  extents_type m_patch_extents{};
  equation_type m_equation{};
};

template <int Rank, typename Localities>
interval_map<octant<Rank>, hpx::id_type>
make_uniform_distribution(int_constant<Rank>, int depth,
                          const Localities& localities) {
  const std::size_t num_locals = localities.size();
  if (num_locals > 0) {
    const std::size_t total_size = index(1) << depth;
    const std::size_t chunk_size = total_size / num_locals;
    std::size_t chunk_rest = total_size % num_locals;
    interval_map<octant<Rank>, hpx::id_type> dist(localities[0]);
    octant<Rank> lower(depth, {0});
    for (const hpx::id_type& locality : localities) {
      try {
        const std::size_t n = chunk_rest ? chunk_size + 1 : chunk_size;
        chunk_rest = chunk_rest - 1;
        const octant<Rank> upper = advance(lower, n);
        dist.insert(lower, upper, locality);
        lower = upper;
      } catch (std::exception&) {
        const octant<Rank> upper = advance(lower, chunk_size - 1);
        dist.insert(lower, upper, locality);
        dist.set_default_value(locality);
      }
    }
    return dist;
  }
  return interval_map<octant<Rank>, hpx::id_type>(hpx::find_here());
}

template <typename Distribution, typename Equation, typename Extents>
grid<Equation, Extents> make_grid(const Distribution& distribution, int depth,
                                  const Extents& extents,
                                  const Equation& equation) {
  static constexpr int rank = Extents::rank;
  grid<Equation, Extents> grid(equation, extents);
  octant<rank> oct(depth, {0});
  const octant<rank> last = octant<rank>().upper_descendant(depth);
  while (oct != last) {
    const hpx::id_type& id = distribution[oct];
    grid.insert(grid.end(), oct, grid_node<Equation, Extents>(id, extents));
    oct = oct.next();
  }
  const hpx::id_type& id = distribution[oct];
  grid.insert(grid.end(), oct, grid_node<Equation, Extents>(id, extents));
  return grid;
}
} // namespace distributed

template <typename Equation, typename Extents>
struct grid_traits<distributed::grid<Equation, Extents>> {
  using grid_type = distributed::grid<Equation, Extents>;
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
  using bind_node = distributed::grid_node<Eq, Ex>;

  using location_type = hpx::id_type;

  template <typename T> using future = hpx::future<T>;

  static constexpr int rank = grid_type::rank;

  template <typename AddressType, AddressType Address, typename T>
  using make_action = hpx::actions::make_action<AddressType, Address, T>;

  template <typename F, typename... Args>
  static decltype(auto) dataflow(F&& f, Args&&... args) {
    return hpx::dataflow(hpx::launch::async, std::forward<F>(f),
                         std::forward<Args>(args)...);
  }

  template <typename F, typename... Args>
  static decltype(auto) dataflow_action(F&& f, Args&&... args) {
    return hpx::async(std::forward<F>(f), std::forward<Args>(args)...);
  }

  /// @brief Returns the octant from a specified value
  static const octant_type& octant(const partition_type& partition) noexcept {
    return partition.first;
  }

  /// @brief Returns the octant for a specified parition
  static const octant_type& octant(partition_type& partition) noexcept {
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

  static hpx::id_type get_location(const partition_type& partition) {
    return partition.second.get_location();
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
    return hpx::when_all(projected).then(
        [initial, binary_op = std::move(binary_op)](hpx::future<std::vector<Projected>> ps) {
          std::vector<Projected> projected_vector = ps.get();
          T value = initial;
          for (Projected& projected : projected_vector) {
            value = fub::invoke(binary_op, value, std::move(projected));
          }
          return value;
        });
  }
};

} // namespace fub

namespace hpx {
namespace serialization {
template <typename Ar, typename Rep, typename Period>
void serialize(Ar& archive, std::chrono::duration<Rep, Period>& t, unsigned) {
  Rep rep = t.count();
  archive & rep;
  t = std::chrono::duration<Rep, Period>{rep};
}
} // namespace serialization
} // namespace hpx

#endif // !FUB_DISTRIBUTED_GRID_HPP
