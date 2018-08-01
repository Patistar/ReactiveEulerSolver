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
#include "fub/p4est/ghost_meta_data.hpp"
#include "fub/p4est/quadrant.hpp"
#include "fub/p4est/tree.hpp"
#include "fub/uniform_cartesian_coordinates.hpp"
#include "fub/variable_data.hpp"
#include "fub/variable_view.hpp"

extern "C" {
#include <mpi.h>
#include <p4est_bits.h>
}

#include <memory>

namespace fub {
inline namespace v1 {
namespace p4est {

struct min_level_t {
  int value;
};
min_level_t min_level(int level) { return min_level_t{level}; }

struct fill_uniform_t {};
static constexpr fill_uniform_t fill_uniform;

/// \ingroup p4est
/// This is a wrapper class around a p4est structure which presents a forest of
/// octrees.
/// \nosubgrouping
template <typename VariableList, int Rank> class grid {
public:
  /// \name Public Types
  /// @{

  using variable_list = VariableList;
  using coordinates_type = uniform_cartesian_coordinates<Rank>;
  using extents_type = dynamic_extents_t<Rank>;
  using patch_data = std::vector<double>;
  using patch_data_view = variable_view<variable_list, double, extents_type>;
  using const_patch_data_view =
      variable_view<variable_list, const double, extents_type>;
  using node_type = quadrant<Rank>;

  /// @}

  /// \name Constructors & Destructors
  /// @{

  /// Constructs an empty grid.
  ///
  /// \throws Nothing.
  grid() = default;

  /// Constructs a new p4est grid.
  ///
  /// This constructor will create only one quadrant at the root level.
  /// When creating the quadrant we will allocate patch data associated to that
  /// quadrant.
  ///
  /// \throws Any exception thrown by the allocator.
  grid(MPI_Comm communicator, connectivity<Rank> conn, variable_list list,
       extents_type extents, coordinates_type coordinates);

  /// @{
  /// Constructs a new p4est structure from a intertree connection info
  /// structure and given intial condition function.
  ///
  /// All quadrants will have a minimum refinement level `min_level`. If
  /// `fill_uniform` is true quadrants will be uniformly distributed accross
  /// processes.
  ///
  /// \throws Any exception thrown by the allocator.
  grid(MPI_Comm communicator, connectivity<Rank> conn, int min_level,
       int fill_uniform, variable_list list, extents_type extents,
       coordinates_type coordinates);

  /// Constructs a new p4est structure from a intertree connection info
  /// structure and given intial condition function.
  ///
  /// All quadrants will have a minimum refinement level `min_level`. If this
  /// constructor is invoked, than `fill_uniform` is always true and quadrants
  /// will be uniformly distributed accross processes.
  ///
  /// \throws Any exception thrown by the allocator.
  grid(MPI_Comm communicator, connectivity<Rank> conn, min_level_t min_level,
       fill_uniform_t fill_uniform, variable_list list, extents_type extents,
       coordinates_type coordinates);
  /// @}

  /// @}

  /// \name Observers
  /// @{

  /// Returns meta data describing all vairables on the grid.
  ///
  /// \throws Nothing.
  variable_list get_variable_list() const noexcept;

  /// Returns the base coordinate function which maps patch indices to euclidean
  /// space coordinates for 0th patch level.
  ///
  /// \throws Nothing.
  coordinates_type get_coordinates() const noexcept;

  /// Returns the coordinate function which maps patch indices to euclidean
  /// space coordinates for a specified quadrant.
  ///
  /// \throws Nothing.
  coordinates_type get_coordinates(const quadrant<Rank>&) const noexcept;

  /// Returns the extents object which every patch is going to have.
  ///
  /// \throws Nothing.
  extents_type get_patch_extents() const noexcept;

  /// Returns a view of all quadrants which are contained in the local forest.
  ///
  /// \return Returns a span to constant nodes.
  span<const tree<Rank>> get_local_trees() const noexcept;

  /// Returns a view of all quadrants which are ghost to the local forest.
  ///
  /// \return Returns a span to constant nodes.
  span<const quadrant<Rank>> get_ghost_quadrants() const noexcept;

  /// Returns a view of all quadrants which are local in this forest and ghost
  /// to some remote forest.
  ///
  /// \return Returns a span to constant nodes.
  span<const quadrant<Rank>> get_mirror_quadrants() const noexcept;

  /// Search for the specified neighbor accross face `face` for the quadrant
  /// `quad`.
  ///
  /// \param[in] quad  the origin quadrant
  /// \param[in] face  the face to look across for a neighbor.
  ///
  /// \return The neighbor quadrant if it is part of the local or ghost mesh.
  /// Otherwise returns an empty optional.
  ///
  /// \throws Nothing.
  optional<quadrant<Rank>> find_face_neighbor(quadrant<Rank> quad,
                                              fub::face face) const noexcept;

  /// Returns true if the speified quadrant `quad` is owned by this process.
  bool is_local(quadrant<Rank> quad) const noexcept;

  /// Returns true if the specified quadrant `quad` is a ghost quadrant.
  bool is_ghost(quadrant<Rank> quad) const noexcept;

  /// Returns true if the specified quadrant `quad` is neither local nor ghost.
  bool is_outside(quadrant<Rank> quad) const noexcept;

  /// Returns true if the specified quadrant `quad` is local for this quadrant
  /// and `quad` is a ghost quadrant for another proccess.
  bool is_mirror(quadrant<Rank> quad) const noexcept;

  /// @}

  /// @{
  /// \name Data Access

  /// @{
  /// Creates a multi-dimensional view onto the patch data which is also
  /// indexable by variable tags.
  ///
  /// If the specified quadrant lies outside of our mesh, we will return an
  /// empty view.
  ///
  /// \param[in] quad  The quadrant which refers to the patch data.
  ///
  /// \return Returns a variable_view which spans all data of a quadrant.
  ///
  /// \post `is_outside(quad) || bool(get_patch_data(quad))`
  ///
  /// \throws Nothing.
  patch_data_view get_patch_data(const quadrant<Rank>& quad) noexcept;

  const_patch_data_view get_patch_data(const quadrant<Rank>& quad) const
      noexcept;
  /// @}

  /// @}

  /// \name Modifiers

  /// Synchronizes the patch data which lie on ghost quadrants.
  ///
  /// This method is a synchronization point and has to be called from every
  /// other process.
  ///
  /// \throws Nothing.
  void exchange_ghost_data();

  /// Refines the local forest.
  ///
  /// This function might allocate new data for the refined quadrants.
  ///
  /// \param[in] predicate Takes a quadrant and signals if a quadrant shall be
  /// refined.
  /// \param[in] transfer This function will be invoked to interpolate data on
  /// the new level.
  ///
  /// \tparam Predicate needs to satisfy `Predicate<quadrant<Rank>>`
  /// \tparam TransferFunction needs to satisfy `Invokable<quadrant<Rank>,
  /// std::array<quadrant<Rank>, 4>>`
  ///
  /// \throws std::bad_alloc on allocation error.
  template <typename Predicate, typename TransferFunction>
  bool refine(Predicate predicate, TransferFunction transfer);

  /// Coarsens the local forest.
  ///
  /// This function might allocate and free data for the coarsened quadrants.
  ///
  /// \param[in] predicate Takes a quadrant and signals if a quadrant shall be
  /// refined.
  /// \param[in] transfer This function will be invoked to interpolate data on
  /// the new level.
  ///
  /// \tparam Predicate needs to satisfy `Predicate<quadrant<Rank>>`
  /// \tparam TransferFunction needs to satisfy `Invokable<quadrant<Rank>,
  /// std::array<quadrant<Rank>, 4>>`
  ///
  /// \throws std::bad_alloc on allocation error.
  template <typename Predicate, typename TransferFunction>
  bool coarsen(Predicate predicate, TransferFunction transfer);

  /// Enforce a 2:1 relationship between neighbors.
  ///
  /// This function might allocate and free data for refined quadrants.
  ///
  /// \param[in] transfer This function will be invoked to interpolate data on
  /// the new level.
  ///
  /// \tparam TransferFunction needs to satisfy `Invokable<quadrant<Rank>,
  /// std::array<quadrant<Rank>, 4>>`
  ///
  /// \throws std::bad_alloc on allocation error.
  template <typename TransferFunction> void balance(TransferFunction transfer);

private:
  variable_list m_list{};
  extents_type m_patch_extents{};
  coordinates_type m_coordinates{};
  std::vector<patch_data> m_patch_datas{};
  std::vector<patch_data> m_ghost_datas{};
  connectivity<Rank> m_connectivity{};
  forest<Rank> m_forest{};
  ghost_meta_data<Rank> m_ghost_meta_data{};
};

template <typename VariableList, int Rank>
grid<VariableList, Rank>::grid(MPI_Comm communicator, connectivity<Rank> conn,
                               variable_list list, extents_type extents,
                               coordinates_type coordinates)
    : grid(communicator, std::move(conn), 0, 1, list, extents, coordinates) {}

template <typename VariableList, int Rank>
grid<VariableList, Rank>::grid(MPI_Comm communicator, connectivity<Rank> conn,
                               int min_level, int fill_uniform,
                               variable_list list, extents_type extents,
                               coordinates_type coordinates)
    : m_list(list), m_patch_extents(extents), m_coordinates(coordinates),
      m_connectivity(std::move(conn)),
      m_forest(communicator, m_connectivity, 0, min_level, fill_uniform,
               [&](auto* g, auto, auto* q) {
                 q->p.user_long = static_cast<long>(m_patch_datas.size());
                 m_patch_datas.emplace_back(m_list.size() *
                                            size(m_patch_extents));
               }),
      m_ghost_meta_data{m_forest} {}

template <typename VariableList, int Rank>
grid<VariableList, Rank>::grid(MPI_Comm communicator, connectivity<Rank> conn,
                               min_level_t min_level, fill_uniform_t,
                               variable_list list, extents_type extents,
                               coordinates_type coordinates)
    : grid(communicator, std::move(conn), min_level.value, 1, list, extents,
           coordinates) {}

template <typename VariableList, int Rank>
typename grid<VariableList, Rank>::variable_list
grid<VariableList, Rank>::get_variable_list() const noexcept {
  return m_list;
}

template <typename VariableList, int Rank>
typename grid<VariableList, Rank>::coordinates_type
grid<VariableList, Rank>::get_coordinates() const noexcept {
  return m_coordinates;
}

template <typename VariableList, int Rank>
typename grid<VariableList, Rank>::coordinates_type
grid<VariableList, Rank>::get_coordinates(const node_type& node) const
    noexcept {
  octant<Rank> oct{node.get_level(), node.get_coordinates()};
  return adapt(m_coordinates, oct);
}

template <typename VariableList, int Rank>
typename grid<VariableList, Rank>::extents_type
grid<VariableList, Rank>::get_patch_extents() const noexcept {
  return m_patch_extents;
}

template <typename VariableList, int Rank>
span<const tree<Rank>> grid<VariableList, Rank>::get_local_trees() const
    noexcept {
  const tree<Rank>* pointer =
      reinterpret_cast<const tree<Rank>*>(m_forest->trees->array);
  std::ptrdiff_t size = m_forest->trees->elem_count;
  return {pointer, size};
}

template <typename VariableList, int Rank>
span<const quadrant<Rank>> grid<VariableList, Rank>::get_ghost_quadrants() const
    noexcept {
  const quadrant<Rank>* pointer =
      reinterpret_cast<const quadrant<Rank>*>(m_ghost_meta_data->ghosts.array);
  std::ptrdiff_t size = m_ghost_meta_data->ghosts.elem_count;
  return {pointer, size};
}

template <typename VariableList, int Rank>
span<const quadrant<Rank>>
grid<VariableList, Rank>::get_mirror_quadrants() const noexcept {
  const quadrant<Rank>* pointer =
      reinterpret_cast<const quadrant<Rank>*>(m_ghost_meta_data->mirrors.array);
  std::ptrdiff_t size = m_ghost_meta_data->mirrors.elem_count;
  return {pointer, size};
}

template <typename VariableList, int Rank>
typename grid<VariableList, Rank>::patch_data_view
grid<VariableList, Rank>::get_patch_data(const node_type& node) noexcept {
  patch_data& data = m_patch_datas[node.index()];
  return patch_data_view(get_variable_list(), data, get_patch_extents());
}

template <typename VariableList, int Rank>
typename grid<VariableList, Rank>::const_patch_data_view
grid<VariableList, Rank>::get_patch_data(const node_type& node) const noexcept {
  const patch_data& data = m_patch_datas[node.index()];
  return const_patch_data_view(get_variable_list(), data, get_patch_extents());
}

template <typename VariableList, int Rank>
optional<quadrant<Rank>>
grid<VariableList, Rank>::find_face_neighbor(quadrant<Rank> quad, face f) const
    noexcept {
  int face_id = as_int(f);
  quadrant<Rank> face_quad = face_neighbor(quad, face_id);
  int owner = -1;
  auto idx =
      face_quadrant_exists(m_forest, m_ghost_meta_data, 0,
                           &face_quad.get_native(), &face_id, nullptr, &owner);
  if (idx >= 0) {
    return face_quad;
  }
  // If we can not find the neighbor we assert that the requested neighbor is on
  // the domain boundary.
  assert(idx == -2);
  return {};
}

template <typename VariableList, int Rank>
void grid<VariableList, Rank>::exchange_ghost_data() {}

} // namespace p4est
} // namespace v1
} // namespace fub

#endif