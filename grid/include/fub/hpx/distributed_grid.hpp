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

#ifndef DISTRIBUTED_GRID_HPP
#define DISTRIBUTED_GRID_HPP

#if defined(__CLANG__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#endif
#include <hpx/include/components.hpp>
#if defined(__CLANG__)
#pragma clang diagnostic pop
#endif

#include "fub/interval_map.hpp"
#include "fub/octree.hpp"
#include "fub/optional.hpp"
#include "fub/patch.hpp"

namespace fub {
namespace hpx {

template <typename Eq, typename E> struct distributed_grid_patch {
  using type = patch<as_tuple_t<typename Eq::complete_state>, E>;
};

template <typename Eq, typename Extents>
class distributed_grid_node : public ::hpx::components::simple_component_base<
                                  distributed_grid_node<Eq, Extents>> {
public:
  using patch_type = typename distributed_grid_patch<Eq, Extents>::type;
  distributed_grid_node() = default;
  explicit distributed_grid_node(patch_type&& patch) noexcept
      : m_patch{std::move(patch)} {}
  explicit distributed_grid_node(const patch_type& patch) noexcept
      : m_patch{patch} {}

  patch_type get_patch() const { return m_patch; }

private:
  patch_type m_patch;
};

template <typename Equation, typename PatchExtents> class distributed_grid {
public:
  static constexpr int rank = PatchExtents::rank;

  using equation_type = Equation;
  using extents_type = PatchExtents;
  using patch_type =
      patch<as_tuple_t<typename Equation::complete_state>, extents_type>;
  using node_type = distributed_grid_node<equation_type, extents_type>;
  using mapped_type = ::hpx::shared_future<node_type>;
  using octree_type = octree<mapped_type, rank>;
  using octant_type = octant<rank>;
  using partition_type = typename octree_type::value_type;
  using iterator = typename octree_type::iterator;
  using const_iterator = typename octree_type::const_iterator;

  distributed_grid() = default;

  // Accessors

  octree_type octree() const noexcept;
  extents_type patch_extents() const noexcept;
  equation_type equation() const noexcept;

  // Visitor

  std::function<void(partition_type)>
  foreach_partition(std::function<void(partition_type)> visitor);

  // Modifiers

  bool insert(const octant<rank>& octant, mapped_type&& value);
  mapped_type find(const octant<rank>&);

  // Iterators

private:
  octree_type m_tree{};
  extents_type m_patch_extents{};
  equation_type m_equation{};
};

template <typename Eq, typename E>
class distributed_grid_client
    : ::hpx::components::client_base<distributed_grid_client<Eq, E>,
                                     distributed_grid<Eq, E>> {};

} // namespace hpx
} // namespace fub

#endif // !DISTRIBUTED_GRID_HPP
