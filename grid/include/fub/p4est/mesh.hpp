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

/// \defgroup p4est
/// This module contains all wrapper types which interface to the p4est library.

#ifndef FUB_P4EST_MESH_HPP
#define FUB_P4EST_MESH_HPP

#include "fub/p4est/forest.hpp"
#include "fub/p4est/ghost_layer.hpp"

#include "fub/span.hpp"

#include <boost/container/pmr/vector.hpp>
#include <boost/container/static_vector.hpp>

namespace fub {
inline namespace v1 {
namespace p4est {

template <int Rank> class mesh;

template <int Rank> struct face_info_data {
  quadrant<Rank> quad;
  bool is_ghost;
};

template <int Rank> class face_info {
public:
  face_info(face_info_data<Rank> coarse,
            std::array<face_info_data<Rank>, 2> fine, fub::face face) noexcept
      : m_coarse{coarse}, m_fine{fine}, m_face{face} {}

  void data(span<byte> new_data) noexcept { m_data = new_data; }
  span<byte> data() noexcept { return m_data; }
  span<const byte> data() const noexcept { return m_data; }

  const face_info_data<Rank>& coarse() const noexcept { return m_coarse; }

  const std::array<face_info_data<Rank>, 2>& fine() const noexcept {
    return m_fine;
  }

  fub::face face() const noexcept { return m_face; }

private:
  face_info_data<Rank> m_coarse;
  std::array<face_info_data<Rank>, 2> m_fine;
  fub::face m_face;
  span<byte> m_data;
};

struct mesh_data_size {
  std::ptrdiff_t local;
  std::ptrdiff_t ghost;
  std::ptrdiff_t coarse_fine_interface;
};

template <> class mesh<2> {
public:
  using allocator_type = boost::container::pmr::polymorphic_allocator<byte>;

  using face_neighbors_t =
      std::array<boost::container::static_vector<quadrant<2>, 4>, 4>;

  using projection_type =
      function_ref<void(const quadrant<2>&, span<const byte>, span<byte>)>;

  mesh(forest<2> forest, ghost_layer<2> ghost, mesh_data_size size,
       allocator_type alloc = allocator_type());

  mesh(const mesh&);
  mesh& operator=(const mesh&);

  mesh(mesh&&) noexcept;
  mesh& operator=(mesh&&) noexcept;

  ~mesh() noexcept;

  void allocate(mesh_data_size size, allocator_type alloc = allocator_type());

  void deallocate() noexcept;

  void reset(forest<2> forest, ghost_layer<2> ghost, mesh_data_size size,
             allocator_type alloc = allocator_type());

  int synchronize_ghost_layer(projection_type projection);

  const forest<2>& forest() const noexcept;

  const ghost_layer<2>& ghost_layer() const noexcept;

  span<byte> local_data(quadrant<2> quad) noexcept;

  span<const byte> local_data(quadrant<2> quad) const noexcept;

  span<byte> local_data(std::ptrdiff_t idx) noexcept;

  span<const byte> local_data(std::ptrdiff_t idx) const noexcept;

  span<byte> ghost_data(quadrant<2> quad) noexcept;

  span<const byte> ghost_data(quadrant<2> quad) const noexcept;

  span<const quadrant<2>> face_neighbors(quadrant<2> quad, face f) const
      noexcept;

  span<const std::ptrdiff_t> quadrants_at_level(int level) const noexcept;

  span<const face_info<2>> coarse_fine_interfaces() const noexcept;

private:
  fub::p4est::forest<2> m_forest;

  fub::p4est::ghost_layer<2> m_ghost;

  allocator_type m_allocator;

  /// For each local quadrant a data buffer.
  span<span<byte>> m_local_to_data;

  /// For each ghost quadrant a data buffer.
  span<span<byte>> m_ghost_to_data;

  /// For each mirror quadrant a data buffer.
  span<span<byte>> m_mirror_to_data;

  /// For each local quadrant its face neighbor relationships.
  span<face_neighbors_t> m_quad_to_face_neighbors;

  /// For each refinement level a list of local quadrant indices.
  std::array<span<std::ptrdiff_t>, P4EST_QMAXLEVEL> m_level_to_local;

  /// A list of "hanging" faces between quadrants at different refinement
  /// levels.
  span<face_info<2>> m_coarse_fine_interfaces;
};

boost::container::static_vector<quadrant<2>, 4>
find_children(const forest<2>& f, const quadrant<2>& q) noexcept;

} // namespace p4est
} // namespace v1
} // namespace fub

#endif