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
#include "fub/p4est/unique_quadrant_data.hpp"
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

template <typename T> struct MPI_datatype;

template <> struct MPI_datatype<float> {
  operator MPI_Datatype() { return MPI_FLOAT; }
};

template <> struct MPI_datatype<double> {
  operator MPI_Datatype() { return MPI_DOUBLE; }
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

  /// The mesh type
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
    coordinates_type coordinates;
  };

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
  span<const mesh_quadrant<Rank>> face_neighbors(quadrant<Rank> quad,
                                                 face f) const noexcept;

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
  template <typename Projection, typename Projected>
  int async_exchange_ghost_data(Projection projection,
                                span<Projected> ghost_data,
                                span<Projected> mirror_data,
                                span<MPI_Request> requests) const;

  int exchange_quadrant_data();

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
  forest_type m_forest;
  ghost_layer_type m_ghost_layer;
  patch_data_info m_patch_data_info;
  mesh_type m_mesh;
  unique_quadrant_data<floating_point[]> m_data;
  unique_quadrant_data<floating_point[]> m_ghost_data;
};

////////////////////////////////////////////////////////////////////////////////
//                                                               Implementation
// CONSTRUCTOR

template <typename VariableList, int Rank, typename FloatingPointType>
basic_grid<VariableList, Rank, FloatingPointType>::basic_grid(
    forest_type forest, ghost_layer_type ghost_layer, patch_data_info info)
    : m_forest{std::move(forest)}, m_ghost_layer{std::move(ghost_layer)},
      m_patch_data_info{info}, m_mesh{m_forest, m_ghost_layer},
      m_data{m_forest.local_num_quadrants(),
             patch_view::static_size(m_patch_data_info.variable_list,
                                     m_patch_data_info.patch_extents)} {}

////////////////////////////////////////////////////////////////////////////////
//                                                Access Patch Data Descriptors

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

////////////////////////////////////////////////////////////////////////////////
//                                                            Access Patch Data

template <typename VariableList, int Rank, typename FloatingPointType>
typename basic_grid<VariableList, Rank, FloatingPointType>::const_patch_view
basic_grid<VariableList, Rank, FloatingPointType>::patch_data(
    int treeidx, int localnum) const noexcept {
  const std::ptrdiff_t index = m_forest.trees()[treeidx].offset() + localnum;
  const span<const floating_point> data = m_data[index];
  return const_patch_view(variable_list(), data, patch_extents());
}

template <typename VariableList, int Rank, typename FloatingPointType>
typename basic_grid<VariableList, Rank, FloatingPointType>::patch_view
basic_grid<VariableList, Rank, FloatingPointType>::patch_data(
    int treeidx, int localnum) noexcept {
  const std::ptrdiff_t index = m_forest.trees()[treeidx].offset() + localnum;
  const span<floating_point> data = m_data[index];
  return patch_view(variable_list(), data, patch_extents());
}

template <typename VariableList, int Rank, typename FloatingPointType>
typename basic_grid<VariableList, Rank, FloatingPointType>::const_patch_view
basic_grid<VariableList, Rank, FloatingPointType>::patch_data(
    quadrant<Rank> quad) const noexcept {
  if (is_local(quad)) {
    return patch_data(quad.which_tree(), quad.local_num());
  }
  return {};
}

template <typename VariableList, int Rank, typename FloatingPointType>
typename basic_grid<VariableList, Rank, FloatingPointType>::patch_view
basic_grid<VariableList, Rank, FloatingPointType>::patch_data(
    quadrant<Rank> quad) noexcept {
  if (is_local(quad)) {
    return patch_data(quad.which_tree(), quad.local_num());
  }
  return {};
}

/// Returns the underlying forest object.
template <typename VariableList, int Rank, typename FloatingPointType>
const typename basic_grid<VariableList, Rank, FloatingPointType>::forest_type&
basic_grid<VariableList, Rank, FloatingPointType>::forest() const noexcept {
  return m_forest;
}

/// Returns the ghost layer object.
template <typename VariableList, int Rank, typename FloatingPointType>
const typename basic_grid<VariableList, Rank,
                          FloatingPointType>::ghost_layer_type&
basic_grid<VariableList, Rank, FloatingPointType>::ghost_layer() const
    noexcept {
  return m_ghost_layer;
}

/// \name Quadrant Observers

/// Returns all quadrants which share a face with the specified quadrant quad.
///
/// \throws Nothing.
template <typename VariableList, int Rank, typename FloatingPointType>
span<const mesh_quadrant<Rank>>
basic_grid<VariableList, Rank, FloatingPointType>::face_neighbors(
    quadrant<Rank> quad, face f) const noexcept {
  std::ptrdiff_t locidx =
      forest().trees()[quad.which_tree()].offset() + quad.local_num();
  return m_mesh.face_neighbors(locidx, f);
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
template <typename Projection, typename Projected>
int basic_grid<VariableList, Rank, FloatingPointType>::
    async_exchange_ghost_data(Projection projection, span<Projected> ghost_data,
                              span<Projected> mirror_data,
                              span<MPI_Request> requests) const {
  assert(mirror_data.size() >= m_ghost_layer.mirrors().size());
  assert(ghost_data.size() >= m_ghost_layer.quadrants().size());
  assert(requests.size() >= mirror_data.size() + ghost_data.size());
  const std::ptrdiff_t mirror_size = m_ghost_layer.mirrors_by_process().size();
  std::ptrdiff_t req_counter = 0;
  std::ptrdiff_t rank = 0;
  const MPI_Comm comm = m_forest.mpi_communicator();
  const span<const quadrant<Rank>> mirrors = m_ghost_layer.mirrors();
  const span<const int> mirrors_by_proc = m_ghost_layer.mirrors_by_process();
  const span<const int> mirrors_by_proc_offsets = m_ghost_layer.mirrors_by_process_offsets();
  for (std::ptrdiff_t i = 0; i < mirror_size; ++i) {
    while (i >= mirrors_by_proc_offsets[rank + 1]) {
      rank += 1;
    }
    assert(rank < m_forest.mpi_size());
    const std::ptrdiff_t mirror_idx = mirrors_by_proc[i];
    const quadrant<Rank> mirror_quad = mirrors[mirror_idx];
    const const_patch_view data = patch_data(mirror_quad);
    fub::invoke(projection, data, mirror_data[mirror_idx]);
    using T = remove_cvref_t<
        std::remove_pointer_t<decltype(mirror_data[mirror_idx].data())>>;
    const int ec =
        MPI_Isend(mirror_data[mirror_idx].data(),
                  mirror_data[mirror_idx].size(), MPI_datatype<T>{}, rank,
                  mirror_quad.local_num(), comm, &requests[req_counter++]);
    if (ec != MPI_SUCCESS) {
      return ec;
    }
  }
  span<const quadrant<Rank>> ghosts = m_ghost_layer.quadrants();
  span<const int> ghosts_offsets = m_ghost_layer.quadrants_process_offsets();
  std::ptrdiff_t ghost_size = ghosts.size();
  rank = 0;
  for (std::ptrdiff_t ghost_idx = 0; ghost_idx < ghost_size; ++ghost_idx) {
    while (ghost_idx >= ghosts_offsets[rank + 1]) {
      rank += 1;
    }
    assert(rank < m_forest.mpi_size());
    using T = remove_cvref_t<
        std::remove_pointer_t<decltype(ghost_data[ghost_idx].data())>>;
    const quadrant<Rank> ghost_quad = ghosts[ghost_idx];
    const int tag = ghost_quad.local_num();
    const int ec =
        MPI_Irecv(ghost_data[ghost_idx].data(), ghost_data[ghost_idx].size(),
                  MPI_datatype<T>{}, rank, tag, comm, &requests[req_counter++]);
    if (ec != MPI_SUCCESS) {
      return ec;
    }
  }
  return MPI_SUCCESS;
}

template <typename VariableList, int Rank, typename FloatingPointType>
int basic_grid<VariableList, Rank,
               FloatingPointType>::exchange_quadrant_data() {
  m_ghost_data = unique_quadrant_data<FloatingPointType[]>(
      m_ghost_layer.quadrants().size(),
      patch_view::static_size(m_patch_data_info.variable_list,
                              m_patch_data_info.patch_extents),
      m_data.get_allocator());
  unique_quadrant_data<FloatingPointType[]> mirror_data(
      m_ghost_layer.mirrors().size(),
      patch_view::static_size(m_patch_data_info.variable_list,
                              m_patch_data_info.patch_extents),
      m_data.get_allocator());
  std::vector<MPI_Request> requests(m_ghost_layer.quadrants().size() +
                                    m_ghost_layer.mirrors_by_process().size());
  async_exchange_ghost_data(
      [](const_patch_view origin, span<FloatingPointType> mirror) {
        std::memcpy(mirror.data(), origin.span().data(),
                    origin.span().byte_size());
      },
      m_ghost_data.data(), mirror_data.data(), requests);
  return MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
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
  fub::p4est::ghost_layer<Rank> new_ghost_layer(new_forest);
  fub::p4est::mesh<Rank> new_mesh(new_forest, new_ghost_layer);
  unique_quadrant_data<floating_point[]> new_data(
      new_forest.local_num_quadrants(),
      patch_view::static_size(m_patch_data_info.variable_list,
                              m_patch_data_info.patch_extents),
      m_data.get_allocator());

  span<const tree<Rank>> trees = forest().trees();
  for (const tree<Rank>& tree : trees) {
    span<const quadrant<Rank>> old_quadrants = tree.quadrants();
    for (quadrant<Rank> quad : old_quadrants) {
      /// Check if quad is not refined.
      optional<std::ptrdiff_t> idx = find(new_forest, quad.which_tree(), quad);
      if (idx) {
        const_patch_view out_view = patch_data(quad);
        span<const quadrant<Rank>> out_quad{&quad, 1};
        span<const_patch_view> out_data{&out_view, 1};
        span<floating_point> data = new_data[*idx];
        patch_view in_view(variable_list(), data, patch_extents());
        span<patch_view> in_data{&in_view, 1};
        const_replace_data<basic_grid> out{out_quad, out_data};
        replace_data<basic_grid> in{out_quad, in_data};
        fub::invoke(replace, out, in);
      }

      /// Check if quad was refined.
      auto children = fub::p4est::children(quad);
      idx = find(new_forest, quad.which_tree(), children[0]);
      if (idx) {
        const_patch_view out_view = patch_data(quad);
        span<const quadrant<Rank>> out_quad{&quad, 1};
        span<const_patch_view> out_data{&out_view, 1};
        boost::container::static_vector<patch_view, 8> new_quad_data{};
        span<floating_point> data = new_data[*idx];
        patch_view view(variable_list(), data, patch_extents());
        new_quad_data.push_back(view);
        for (int i = 1; i < children.size(); ++i) {
          idx = find(new_forest, quad.which_tree(), children[i]);
          assert(idx);
          data = new_data[*idx];
          patch_view view(variable_list(), data, patch_extents());
          new_quad_data.push_back(view);
        }
        const_replace_data<basic_grid> out{out_quad, out_data};
        replace_data<basic_grid> in{children, new_quad_data};
        fub::invoke(replace, out, in);
      }
    }
  }
  m_forest = std::move(new_forest);
  m_ghost_layer = std::move(new_ghost_layer);
  m_mesh = std::move(new_mesh);
  m_data = std::move(new_data);
}

} // namespace p4est
} // namespace v1
} // namespace fub

#endif