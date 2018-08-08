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

#ifndef FUB_P4EST_GRID_HPP
#define FUB_P4EST_GRID_HPP

#include "fub/functional.hpp"
#include "fub/p4est/connectivity.hpp"
#include "fub/p4est/forest.hpp"
#include "fub/p4est/ghost_layer.hpp"
#include "fub/p4est/mesh.hpp"
#include "fub/p4est/quadrant.hpp"
#include "fub/p4est/tree.hpp"
#include "fub/uniform_cartesian_coordinates.hpp"
#include "fub/variable_data.hpp"
#include "fub/variable_view.hpp"

extern "C" {
#include <mpi.h>
#include <p4est_bits.h>
}

#include <array>
#include <cstring>
#include <map>
#include <memory>

namespace fub {
inline namespace v1 {
namespace p4est {

template <typename Grid> struct replace_data {
  span<const typename Grid::quadrant_type> quadrants;
  span<typename Grid::patch_view> datas;
};

template <typename Grid> struct const_replace_data {
  span<const typename Grid::quadrant_type> quadrants;
  span<typename Grid::const_patch_view> datas;
};

/// \ingroup p4est
/// This class holds a patch data with regard to a fixed forest and ghost
/// strucutre.
template <typename VariableList, int Rank, typename FloatingPointType>
class basic_grid {
public:
  /// \name Public Types

  /// The value type is an arithmetic type which is currently only double or
  /// float.
  using floating_point = FloatingPointType;

  static_assert(std::is_floating_point<floating_point>::value,
                "FloatingPointType must be a floating point type.");

  /// The variable list provide named access to data patches.
  using variable_list_type = VariableList;

  static_assert(is_variable_list<VariableList>::value,
                "VariableList is not a variable list type.");

  /// The coordinate mapping which embed patches into the cartesian space.
  using coordinates_type = uniform_cartesian_coordinates<Rank>;

  /// The uniform extents for each single patch.
  using extents_type = dynamic_extents_t<Rank>;

  /// The p4est quadrant type.
  using quadrant_type = quadrant<Rank>;

  /// The p4est forest type.
  using forest_type = forest<Rank>;

  /// The p4est ghost layer type.
  using ghost_layer_type = ghost_layer<Rank>;

  using mesh_type = mesh<Rank>;

  /// A view to mutable patch data.
  using patch_view =
      variable_view<variable_list_type, floating_point, extents_type>;

  /// A view to const patch data.
  using const_patch_view =
      variable_view<variable_list_type, const floating_point, extents_type>;

  struct patch_data_info {
    variable_list_type variable_list;
    extents_type patch_extents;
    extents_type ghost_extents;
    coordinates_type coordinates;
  };

  static std::ptrdiff_t local_data_byte_size(patch_data_info info) noexcept;
  static std::ptrdiff_t ghost_data_byte_size(patch_data_info info) noexcept;

  /// \name Constructors, Assignment & Destructors

  /// Constructs a new grid with zero initialised patch datas.
  ///
  /// \param[in] forest  A forest object which describes how many data patches
  /// we are going to allocate.
  ///
  /// \param[in] list  A variable_list object to describe the variables on the
  /// patch data.
  ///
  /// \param[in] extents  The index extents for each patch.
  ///
  /// \param[in] coordinates  A coordinate mapping which maps the level 0 patch
  /// to the whole domain.
  ///
  /// \throws std::bad_alloc  On Allocation error.
  basic_grid(forest_type forest, ghost_layer_type ghost,
             patch_data_info patch_data_info);

  /// \name Patch Data Observers

  /// Returns meta data describing all vairables on the grid.
  ///
  /// \throws Nothing.
  variable_list_type variable_list() const noexcept;

  /// Returns the base coordinate function which maps patch indices to euclidean
  /// space coordinates for 0th patch level.
  ///
  /// \throws Nothing.
  coordinates_type coordinates() const noexcept;

  /// Returns the coordinate function which maps patch indices to euclidean
  /// space coordinates for a specified quadrant.
  ///
  /// \throws Nothing.
  coordinates_type coordinates(quadrant<Rank>) const noexcept;

  /// Returns the extents object which every patch is going to have.
  ///
  /// \throws Nothing.
  extents_type patch_extents() const noexcept;

  /// Returns the extents object which every ghost is going to have.
  ///
  /// \throws Nothing.
  extents_type ghost_extents() const noexcept;

  /// \name Patch Data Access

  /// @{
  /// Creates a multi-dimensional view of the patch data.
  ///
  /// If the specified quadrant lies outside of our mesh, the function will
  /// return an empty view.
  ///
  /// \param[in] quad  The quadrant which refers to the patch data.
  ///
  /// \return Returns a variable_view which spans all data of a quadrant.
  ///
  /// \post `is_outside(quad) || bool(patch_data(quad))`
  ///
  /// \throws Nothing.
  const_patch_view patch_data(quadrant<Rank> quad) const noexcept;
  patch_view patch_data(quadrant<Rank> quad) noexcept;
  /// @}

  /// @{
  /// Creates a multi-dimensional view of the patch data.
  ///
  /// If the specified quadrant lies outside of our mesh, the function will
  /// return an empty view.
  ///
  /// \param[in] which_tree  The tree index which the quadrant belongs to.
  ///
  /// \param[in] local_idx  The local quadrant index in the tree.
  ///
  /// \return Returns a variable_view which spans all data of a quadrant.
  ///
  /// \post `is_outside(quad) || bool(patch_data(quad))`
  ///
  /// \throws Nothing.
  const_patch_view patch_data(int which_tree, int local_idx) const noexcept;
  patch_view patch_data(int which_tree, int local_idx) noexcept;
  /// @}

  /// \name p4est related Member Access

  /// Returns the underlying forest object.
  const forest_type& forest() const noexcept;

  /// Returns the ghost layer object.
  const ghost_layer_type& ghost_layer() const noexcept;

  /// \name Quadrant Observers

  /// Returns all quadrants which share a face with the specified quadrant quad.
  ///
  /// \throws Nothing.
  span<const quadrant<Rank>> face_neighbors(quadrant<Rank> quad, face f) const
      noexcept;

  /// Returns the owner rank of the specfied ghost quadrant
  optional<int> find_owner_mpi_rank(quadrant<Rank> quad) const noexcept;

  /// Returns true if the speified quadrant `quad` is owned by this process.
  bool is_local(quadrant<Rank> quad) const noexcept;

  /// Returns true if the specified quadrant `quad` is a ghost quadrant.
  bool is_ghost(quadrant<Rank> quad) const noexcept;

  /// Returns true if the specified quadrant `quad` is neither local nor ghost.
  bool is_outside(quadrant<Rank> quad) const noexcept;

  /// Returns true if the specified quadrant `quad` is local for this quadrant
  /// and `quad` is a ghost quadrant for another proccess.
  bool is_mirror(quadrant<Rank> quad) const noexcept;

  /// \name Modifiers

  /// Synchronizes the patch data which lie on ghost quadrants for the specified
  /// variables.
  ///
  /// This method is a synchronization point and has to be called from every
  /// other process.
  ///
  /// \param[in] ghost_width  The amount of cells in the `x` direction which get
  /// synchronized.
  ///
  /// \param[in] variables  A list of variables which data shall be send and
  /// recieved. This can be empty to simply snychronize everything.
  ///
  /// \return Returns MPI_SUCCESS on success otherwise an MPI error code.
  ///
  /// \throws std::bad_alloc on allocation failure.
  template <typename... Variables>
  int synchronize_ghost_data(Variables... variables);

  /// Changes the quadrant distribution according to `new_forest`.
  ///
  /// This function invokes the `replace` function whenever an old quadrant was
  /// either refined and replaced by new quadrants, or if a family of quadrants
  /// were coarsened and replaced by one quadrant.
  ///
  /// \param[in] new_forest  describes the new forest of quadrants.
  ///
  /// \param[in] replace  An invokable object which replaces between the
  /// outgoing with incoming quadrants data values.
  ///
  /// \tparam Replace  A function type which satisfies `std::Invokable<void,
  /// span<const_patch_view>, span<patch_view>>`.
  ///
  /// \throws std::bad_alloc on data allocation failure.
  template <typename Replace>
  void transfer_data_to(forest_type new_forest, Replace replace);

private:
  mesh_type m_mesh;
  patch_data_info m_patch_data_info;
};

////////////////////////////////////////////////////////////////////////////////
//                                                               Implementation

template <typename VariableList, int Rank, typename FloatingPointType>
std::ptrdiff_t
basic_grid<VariableList, Rank, FloatingPointType>::local_data_byte_size(
    patch_data_info info) noexcept {
  return std::ptrdiff_t(sizeof(floating_point)) *
         patch_view::static_size(info.variable_list, info.patch_extents);
}

template <typename VariableList, int Rank, typename FloatingPointType>
std::ptrdiff_t
basic_grid<VariableList, Rank, FloatingPointType>::ghost_data_byte_size(
    patch_data_info info) noexcept {
  return std::ptrdiff_t(sizeof(floating_point)) *
         patch_view::static_size(info.variable_list, info.ghost_extents);
}

// CONSTRUCTOR

template <typename VariableList, int Rank, typename FloatingPointType>
basic_grid<VariableList, Rank, FloatingPointType>::basic_grid(
    forest_type forest, ghost_layer_type ghost, patch_data_info info)
    : m_mesh{std::move(forest), std::move(ghost),
             mesh_data_size{local_data_byte_size(info),
                            ghost_data_byte_size(info)}},
      m_patch_data_info{info} {}

////////////////////////////////////////////////////////////////////////////////
//                                                Access Patch Data
//                                                Descriptors

template <typename VariableList, int Rank, typename FloatingPointType>
typename basic_grid<VariableList, Rank, FloatingPointType>::variable_list_type
basic_grid<VariableList, Rank, FloatingPointType>::variable_list() const
    noexcept {
  return m_patch_data_info.variable_list;
}

template <typename VariableList, int Rank, typename FloatingPointType>
typename basic_grid<VariableList, Rank, FloatingPointType>::coordinates_type
basic_grid<VariableList, Rank, FloatingPointType>::coordinates() const
    noexcept {
  return m_patch_data_info.coordinates;
}

template <typename VariableList, int Rank, typename FloatingPointType>
typename basic_grid<VariableList, Rank, FloatingPointType>::coordinates_type
basic_grid<VariableList, Rank, FloatingPointType>::coordinates(
    quadrant<Rank> quad) const noexcept {
  return adapt(m_patch_data_info.coordinates, quad);
}

template <typename VariableList, int Rank, typename FloatingPointType>
typename basic_grid<VariableList, Rank, FloatingPointType>::extents_type
basic_grid<VariableList, Rank, FloatingPointType>::patch_extents() const
    noexcept {
  return m_patch_data_info.patch_extents;
}

template <typename VariableList, int Rank, typename FloatingPointType>
typename basic_grid<VariableList, Rank, FloatingPointType>::extents_type
basic_grid<VariableList, Rank, FloatingPointType>::ghost_extents() const
    noexcept {
  return m_patch_data_info.ghost_extents;
}

////////////////////////////////////////////////////////////////////////////////
//                                                            Access Patch
//                                                            Data

template <typename VariableList, int Rank, typename FloatingPointType>
typename basic_grid<VariableList, Rank, FloatingPointType>::const_patch_view
basic_grid<VariableList, Rank, FloatingPointType>::patch_data(
    quadrant<Rank> quad) const noexcept {
  span<const byte> raw_data{};
  if (is_local(quad)) {
    raw_data = m_mesh.local_data(quad);
  } else if (is_ghost(quad)) {
    raw_data = m_mesh.ghost_data(quad);
  }
  const floating_point* pointer =
      reinterpret_cast<const floating_point*>(raw_data.data());
  std::ptrdiff_t size = raw_data.size() / sizeof(floating_point);
  span<const floating_point> data{pointer, size};
  return const_patch_view(variable_list(), data, patch_extents());
}

template <typename VariableList, int Rank, typename FloatingPointType>
typename basic_grid<VariableList, Rank, FloatingPointType>::patch_view
basic_grid<VariableList, Rank, FloatingPointType>::patch_data(
    quadrant<Rank> quad) noexcept {
  span<byte> raw_data{};
  if (is_local(quad)) {
    raw_data = m_mesh.local_data(quad);
    floating_point* pointer =
        reinterpret_cast<floating_point*>(raw_data.data());
    std::ptrdiff_t size = raw_data.size() / sizeof(floating_point);
    span<floating_point> data{pointer, size};
    return patch_view(variable_list(), data, patch_extents());
  } else if (is_ghost(quad)) {
    raw_data = m_mesh.ghost_data(quad);
    floating_point* pointer =
        reinterpret_cast<floating_point*>(raw_data.data());
    std::ptrdiff_t size = raw_data.size() / sizeof(floating_point);
    span<floating_point> data{pointer, size};
    return patch_view(variable_list(), data, patch_extents());
  }
  return {};
}

/// Returns the underlying forest object.
template <typename VariableList, int Rank, typename FloatingPointType>
const typename basic_grid<VariableList, Rank, FloatingPointType>::forest_type&
basic_grid<VariableList, Rank, FloatingPointType>::forest() const noexcept {
  return m_mesh.forest();
}

/// Returns the ghost layer object.
template <typename VariableList, int Rank, typename FloatingPointType>
const typename basic_grid<VariableList, Rank,
                          FloatingPointType>::ghost_layer_type&
basic_grid<VariableList, Rank, FloatingPointType>::ghost_layer() const
    noexcept {
  return m_mesh.ghost_layer();
}

/// \name Quadrant Observers

/// Returns all quadrants which share a face with the specified quadrant quad.
///
/// \throws Nothing.
template <typename VariableList, int Rank, typename FloatingPointType>
span<const quadrant<Rank>>
basic_grid<VariableList, Rank, FloatingPointType>::face_neighbors(
    quadrant<Rank> quad, face f) const noexcept {
  return m_mesh.face_neighbors(quad, f);
}

/// Returns the owner rank of the specfied ghost quadrant
template <typename VariableList, int Rank, typename FloatingPointType>
optional<int>
basic_grid<VariableList, Rank, FloatingPointType>::find_owner_mpi_rank(
    quadrant<Rank> quad) const noexcept {
  if (is_local(quad)) {
    return forest().mpi_rank();
  }
  if (is_ghost(quad)) {
    const auto offsets = ghost_layer().quadrants_process_offsets();
    const auto proc =
        std::lower_bound(offsets.begin(), offsets.end(), quad.local_num());
    return std::distance(offsets.begin(), proc);
  }
  return {};
}

/// Returns true if the speified quadrant `quad` is owned by this process.
template <typename VariableList, int Rank, typename FloatingPointType>
bool basic_grid<VariableList, Rank, FloatingPointType>::is_local(
    quadrant<Rank> quad) const noexcept {
  const int treeidx = quad.which_tree();
  if (0 <= treeidx && treeidx < forest().trees().size()) {
    const tree<Rank>& tree = forest().trees()[treeidx];
    const int local_num = quad.local_num();
    if (0 <= local_num && local_num < tree.quadrants().size()) {
      return quad == tree.quadrants()[local_num];
    }
  }
  return false;
}

/// Returns true if the specified quadrant `quad` is a ghost quadrant.
template <typename VariableList, int Rank, typename FloatingPointType>
bool basic_grid<VariableList, Rank, FloatingPointType>::is_ghost(
    quadrant<Rank> quad) const noexcept {
  const int ghostidx = quad.local_num();
  if (0 <= ghostidx && ghostidx < ghost_layer().quadrants().size()) {
    return ghost_layer().quadrants()[ghostidx] == quad;
  }
  return false;
}

/// Returns true if the specified quadrant `quad` is neither local nor ghost.
template <typename VariableList, int Rank, typename FloatingPointType>
bool basic_grid<VariableList, Rank, FloatingPointType>::is_outside(
    quadrant<Rank> quad) const noexcept {
  return !is_local(quad) && !is_ghost(quad);
}

/// Returns true if the specified quadrant `quad` is local for this quadrant
/// and `quad` is a ghost quadrant for another proccess.
template <typename VariableList, int Rank, typename FloatingPointType>
bool basic_grid<VariableList, Rank, FloatingPointType>::is_mirror(
    quadrant<Rank> quad) const noexcept {
  const int mirroridx = quad.local_num();
  if (0 <= mirroridx && mirroridx < ghost_layer().mirrors().size()) {
    return ghost_layer().mirrors()[mirroridx] == quad;
  }
  return false;
}

/// \name Modifiers

/// Synchronizes the patch data which lie on ghost quadrants for the specified
/// variables.
///
/// This method is a synchronization point and has to be called from every
/// other process.
///
/// \param[in] ghost_width  The amount of cells in the `x` direction which get
/// synchronized.
///
/// \param[in] variables  A list of variables which data shall be send and
/// recieved. This can be empty to simply snychronize everything.
///
/// \return Returns MPI_SUCCESS on success otherwise an MPI error code.
///
/// \throws std::bad_alloc on allocation failure.
template <typename VariableList, int Rank, typename FloatingPointType>
template <typename... Variables>
int basic_grid<VariableList, Rank, FloatingPointType>::synchronize_ghost_data(
    Variables... variables) {
  return m_mesh.synchronize_ghost_layer([&](const quadrant<2>& mirror,
                                            span<const byte> local,
                                            span<byte> remote) {
    assert(local.size() == remote.size());
    std::memcpy(remote.data(), local.data(), remote.size());
  });
}

/// Changes the quadrant distribution according to `new_forest`.
///
/// This function invokes the `replace` function whenever an old quadrant was
/// either refined and replaced by new quadrants, or if a family of quadrants
/// were coarsened and replaced by one quadrant.
///
/// \param[in] new_forest  describes the new forest of quadrants.
///
/// \param[in] replace  An invokable object which replaces between the
/// outgoing with incoming quadrants data values.
///
/// \tparam Replace  A function type which satisfies `std::Invokable<void,
/// span<const_patch_info_t>, span<patch_info_t>>`.
///
/// \throws std::bad_alloc on data allocation failure.
template <typename VariableList, int Rank, typename FloatingPointType>
template <typename Replace>
void basic_grid<VariableList, Rank, FloatingPointType>::transfer_data_to(
    forest_type new_forest, Replace replace) {
  balance(new_forest);
  fub::p4est::ghost_layer<Rank> ghost_layer(new_forest);
  mesh_data_size size{local_data_byte_size(m_patch_data_info),
                      ghost_data_byte_size(m_patch_data_info)};
  fub::p4est::mesh<Rank> new_mesh(std::move(new_forest), std::move(ghost_layer),
                                  size);
  span<const tree<Rank>> trees = forest().trees();
  for (const tree<Rank>& tree : trees) {
    span<const quadrant<Rank>> old_quadrants = tree.quadrants();
    for (quadrant<Rank> quad : old_quadrants) {
      /// Check if quad is not refined.
      optional<std::ptrdiff_t> idx = find(new_mesh.forest(), quad);
      if (idx) {
        const_patch_view out_view = patch_data(quad);
        span<const quadrant<Rank>> out_quad{&quad, 1};
        span<const_patch_view> out_data{&out_view, 1};
        span<floating_point> data =
            span_cast<floating_point>(new_mesh.local_data(*idx));
        patch_view in_view(variable_list(), data, patch_extents());
        span<patch_view> in_data{&in_view, 1};
        const_replace_data<basic_grid> out{out_quad, out_data};
        replace_data<basic_grid> in{out_quad, in_data};
        fub::invoke(replace, out, in);
      }

      /// Check if quad was refined.
      auto children = fub::p4est::children(quad);
      idx = find(new_mesh.forest(), children[0]);
      if (idx) {
        const_patch_view out_view = patch_data(quad);
        span<const quadrant<Rank>> out_quad{&quad, 1};
        span<const_patch_view> out_data{&out_view, 1};
        boost::container::static_vector<patch_view, 8> new_data{};
        span<floating_point> data =
            span_cast<floating_point>(new_mesh.local_data(*idx));
        patch_view view(variable_list(), data, patch_extents());
        new_data.push_back(view);
        for (int i = 1; i < children.size(); ++i) {
          idx = find(new_mesh.forest(), children[i]);
          assert(idx);
          span<floating_point> data =
              span_cast<floating_point>(new_mesh.local_data(*idx));
          patch_view view(variable_list(), data, patch_extents());
          new_data.push_back(view);
        }
        const_replace_data<basic_grid> out{out_quad, out_data};
        replace_data<basic_grid> in{children, new_data};
        fub::invoke(replace, out, in);
      }
    }
  }
  m_mesh = std::move(new_mesh);
}

} // namespace p4est
} // namespace v1
} // namespace fub

#endif