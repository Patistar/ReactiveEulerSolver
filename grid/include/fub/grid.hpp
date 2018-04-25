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

#ifndef FUB_GRID_HPP
#define FUB_GRID_HPP

#include "fub/equation.hpp"
#include "fub/grid/config.hpp"
#include "fub/octree.hpp"
#include "fub/patch.hpp"
#include "fub/patch_view.hpp"

#if defined(FUB_WITH_POLYMORPHIC_ALLOCATOR)
#include <boost/container/pmr/polymorphic_allocator.hpp>
#endif

#include <memory>

namespace fub {

template <typename G> struct grid_traits {
  using grid_type = G;
  using partition_type = typename G::partition_type;
  using extents_type = typename G::extents_type;
  using octant_type = typename G::octant_type;
  using mapped_type = typename G::mapped_type;
  using patch_type = typename G::patch_type;
  using patch_view_type = decltype(make_view(std::declval<patch_type&>()));
  using node_type = typename G::node_type;

  static constexpr int rank = G::rank;

  /// @brief This function attempts to transform its argument to a patch or just
  /// forwards it otherwise.
  ///
  /// @detail The following transformations will be done
  ///
  ///       partition_type -> patch_type
  ///          mapped_type -> patch_type
  ///            node_type -> patch_type
  ///           node_type* -> patch_type
  ///          patch_type* -> patch_type
  /// otherwise:    T -> T
  // {{{

  // partition_type -> patch_type
  static const patch_type& filter_patch(const partition_type& obj) noexcept {
    return filter_patch(std::get<1>(obj));
  }
  static patch_type& filter_patch(partition_type& obj) noexcept {
    return filter_patch(std::get<1>(obj));
  }
  static patch_type filter_patch(partition_type&& obj) noexcept {
    return filter_patch(std::get<1>(std::move(obj)));
  }

  // mapped_type -> patch_type
  static const patch_type& filter_patch(const mapped_type& obj) noexcept {
    return filter_patch(obj.get());
  }
  static patch_type& filter_patch(mapped_type& obj) noexcept {
    return filter_patch(obj.get());
  }
  static patch_type filter_patch(mapped_type&& obj) noexcept {
    return filter_patch(std::move(obj).get());
  }

  // node_type* -> patch_type
  static const patch_type& filter_patch(const node_type* obj) noexcept {
    return obj->patch;
  }
  static patch_type& filter_patch(node_type* obj) noexcept {
    return obj->patch;
  }

  // node_type -> patch_type
  static const patch_type& filter_patch(const node_type& obj) noexcept {
    return obj.patch;
  }
  static patch_type& filter_patch(node_type& obj) noexcept { return obj.patch; }
  static patch_type filter_patch(node_type&& obj) noexcept {
    return std::move(obj.patch);
  }

  // patch_type* -> patch_type
  static patch_type& filter_patch(patch_type* p) noexcept { return *p; }
  static const patch_type& filter_patch(const patch_type* p) noexcept {
    return *p;
  }

  // Forward Rest
  template <typename T> static const T& filter_patch(const T& obj) noexcept {
    return obj;
  }
  // }}}

  template <typename F, typename... Args>
  static decltype(auto) dataflow(F f, Args&&... args) {
    return fub::invoke(std::move(f), filter_patch(std::forward<Args>(args))...);
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
  static T reduce(const G& grid, T initial, BinaryOp binary_op, Proj proj) {
    for (const partition_type& partition : grid) {
      initial = invoke(binary_op, initial, invoke(proj, partition));
    }
    return initial;
  }
}; // namespace fub

template <typename Eq, typename Extents> struct standard_grid_node {
  using patch_type = patch<as_tuple_t<complete_state_t<Eq>>, Extents>;
  standard_grid_node() = default;
  explicit standard_grid_node(patch_type&& b) noexcept : patch{std::move(b)} {}
  explicit standard_grid_node(const patch_type& b) noexcept : patch{b} {}
  patch_type patch;
};

template <typename Equation, typename Extents> class standard_grid {
public:
  static constexpr int rank = Extents::rank;

  using equation_type = Equation;
  using extents_type = Extents;
  using patch_type = patch<as_tuple_t<complete_state_t<Equation>>, Extents>;
  using node_type = standard_grid_node<Equation, Extents>;
#if defined(FUB_WITH_POLYMORPHIC_ALLOCATOR)
  using allocator_type =
      boost::container::pmr::polymorphic_allocator<node_type>;
#else
  using allocator_type = std::allocator<node_type>;
#endif
  using octant_type = fub::octant<rank>;
  using mapped_type = std::shared_ptr<node_type>;
  using octree_type = fub::octree<mapped_type, rank>;
  using partition_type = typename octree_type::value_type;
  using iterator = typename octree_type::iterator;
  using const_iterator = typename octree_type::const_iterator;

private:
  octree_type m_tree{};
  extents_type m_patch_extents{};
  equation_type m_equation{};

public:
  // Constructors

  standard_grid() = default;

  standard_grid(int depth, const extents_type& extents = extents_type(),
                const equation_type& equation = equation_type())
      : m_tree{}, m_patch_extents{extents}, m_equation{equation} {
    octant_type octant(depth, {{0}});
    octant_type last = octant_type::root().upper_descendant(depth);
    while (octant != last) {
      auto pointer = std::allocate_shared<node_type>(allocator_type());
      m_tree.insert(m_tree.end(), {octant, std::move(pointer)});
      octant = octant.next();
    }
    auto pointer = std::allocate_shared<node_type>(allocator_type());
    m_tree.insert(m_tree.end(), {octant, std::move(pointer)});
  }

  // Iterators

  iterator begin() noexcept { return m_tree.begin(); }
  iterator end() noexcept { return m_tree.end(); }
  const_iterator begin() const noexcept { return m_tree.begin(); }
  const_iterator end() const noexcept { return m_tree.end(); }

  // Access Partition

  /// @brief Returns the partition for a specified octant if found. Otherwise it
  /// returns an empty optional.
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

  iterator insert(const_iterator hint, const octant_type& octant,
                  patch_type&& value) {
    auto pointer =
        std::allocate_shared<node_type>(allocator_type(), std::move(value));
    return m_tree.insert(hint, std::make_pair(octant, std::move(pointer)));
  }

  iterator insert(const_iterator hint, const octant_type& octant,
                  const patch_type& value) {
    auto pointer = std::allocate_shared<node_type>(allocator_type(), value);
    return m_tree.insert(hint, std::make_pair(octant, std::move(pointer)));
  }

  std::pair<iterator, bool> insert(const octant_type& octant,
                                   patch_type&& value) {
    auto pointer =
        std::allocate_shared<node_type>(allocator_type(), std::move(value));
    return m_tree.insert(std::make_pair(octant, std::move(pointer)));
  }

  std::pair<iterator, bool> insert(const octant_type& octant,
                                   const patch_type& value) {
    auto pointer = std::allocate_shared<node_type>(allocator_type(), value);
    return m_tree.insert(std::make_pair(octant, std::move(pointer)));
  }

  template <typename... Args>
  std::pair<iterator, bool> emplace(const octant_type& octant, Args&&... args) {
    auto pointer = std::allocate_shared<node_type>(
        allocator_type(), DataPatch(std::forward<Args>(args)...));
    return m_tree.insert(std::make_pair(octant, std::move(pointer)));
  }
};

} // namespace fub

#endif // !GRID_HPP
