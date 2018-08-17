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
#include "fub/face.hpp"

#include "fub/span.hpp"

#include <boost/container/pmr/vector.hpp>
#include <boost/container/static_vector.hpp>

#include <vector>

namespace fub {
inline namespace v1 {
namespace p4est {

/// This class manages objects of type T for a forest and its ghost layer.
template <int Rank> class mesh;

template <int Rank> struct mesh_quadrant {
  quadrant<Rank> quad;
  bool is_ghost;
};

template <typename Q> struct face_and {
  fub::face face;
  Q quadrant;
};

/// A face_info object holds a data object which is mapped to a specified
/// coarse-fine interface.
template <int Rank> struct coarse_fine_interface {
  face_and<mesh_quadrant<Rank>> coarse;
  face_and<std::array<mesh_quadrant<Rank>, 2>> fine;
};

template <> class mesh<2> {
public:
  using face_neighbors_t =
      std::array<boost::container::static_vector<mesh_quadrant<2>, 4>, 4>;

  mesh() = default;

  /// Constructs a mesh from a given forest and ghost layer.
  ///
  /// \note The forest and ghost layer should outlive the mesh.
  mesh(forest<2>& forest, ghost_layer<2>& ghost_layer);

  static constexpr int rank() noexcept { return 2; }

  /// Returns a either up to two quadrants which are face neighbors to `quad`
  /// over face `f`.
  span<const mesh_quadrant<2>> face_neighbors(std::ptrdiff_t index, face f) const
      noexcept {
    return m_quad_to_face_neighbors[index][f];
  }

  /// Returns a list of all coarse fine interface information objects.
  span<const coarse_fine_interface<2>> coarse_fine_interfaces() const
      noexcept {
    return m_coarse_fine_interfaces;
  }

private:
  /// For each local quadrant its face neighbor relationships.
  std::vector<face_neighbors_t> m_quad_to_face_neighbors;

  /// A vector of coarse fine interfaces, sorted by local index.
  std::vector<coarse_fine_interface<2>> m_coarse_fine_interfaces;
};

} // namespace p4est
} // namespace v1
} // namespace fub

#endif