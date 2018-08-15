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

#ifndef FUB_P4EST_FOREST_HPP
#define FUB_P4EST_FOREST_HPP

#include "fub/functional.hpp"
#include "fub/optional.hpp"
#include "fub/span.hpp"

#include "fub/p4est/connectivity.hpp"
#include "fub/p4est/quadrant.hpp"
#include "fub/p4est/tree.hpp"

extern "C" {
#include <p4est.h>
#include <p4est_extended.h>
#include <p8est.h>
}

#include <memory>

namespace fub {
inline namespace v1 {
namespace p4est {

template <int Rank> class forest;

/// \ingroup p4est
/// This is a wrapper for the 2-dimensional forest type `p4est_t`.
template <> class forest<2> {
public:
  /// \name Constructors & Assignment

  /// Constructs an invalid forest.
  forest() = default;

  /// Copy constructor makes a deep copy of the forest but does not copy
  /// user-data in the quadrants.
  ///
  /// \param[in] other  The forest which will be copied.
  ///
  /// \throws std::bad_alloc on allocation error.
  forest(const forest& other);

  /// Copy Assignment makes a deep copy of the forest but does not copy
  /// user-data in the quadrants.
  ///
  /// \param[in] other  The forest which will be copied.
  ///
  /// \throws std::bad_alloc on allocation error.
  forest& operator=(const forest& other);

  /// Move constructor leaves an invalid forest.
  ///
  /// \throws Nothing.
  forest(forest&&) noexcept = default;

  /// Move Assignment leaves an invalid forest.
  ///
  /// \throws Nothing.
  forest& operator=(forest&&) noexcept = default;

  /// Take ownership of the specified p4est pointer.
  ///
  /// \throws Nothing.
  explicit forest(p4est_t* pointer) noexcept;

  /// Constructs a new forest.
  ///
  /// The new forest consists of equi-partitioned root quadrants. When there are
  /// more processors than trees, some processors are empty.
  ///
  /// \param[in] communicator  A MPI communicator which refers to the group of
  /// procecsses to share trees in the forest with.
  ///
  /// \param[in] conn  A connectivity object which describes tree connections in
  /// the forest.
  ///
  /// \throws Nothing.
  ///
  /// \note The connectivity structure must not be destroyed during the lifetime
  /// of this forest.
  forest(MPI_Comm communicator, const connectivity<2>& conn) noexcept;

  /// Constructs a new forest.
  ///
  /// \param[in] communicator  A MPI communicator which refers to the group of
  /// procecsses to share trees in the forest with.
  ///
  /// \param[in] conn  A connectivity object which describes tree connections in
  /// the forest.
  ///
  /// \param[in] min_quads  Minimum initial quadrants per processor. Makes the
  /// refinement pattern mpisize-specific.
  ///
  /// \param[in] min_level  The forest is refined at least to this level. May be
  /// negative or 0, then it has no effect.
  ///
  /// \param[in] fill_uniform  If true, fill the forest with a uniform mesh
  /// instead of the coarsest possible one. The latter is partition-specific so
  /// that is usually not a good idea.
  ///
  /// \throws Nothing.
  ///
  /// \note The connectivity structure must not be destroyed during the lifetime
  /// of this forest.
  forest(MPI_Comm communicator, ::fub::p4est::connectivity<2>& conn,
         int min_quads, int min_level, int fill_uniform) noexcept;

  /// \name Member Accessors

  /// Returns the MPI communicator.
  MPI_Comm mpi_communicator() const noexcept;

  /// Returns this process's MPI rank.
  int mpi_rank() const noexcept;

  /// Returns the number of MPI processes.
  int mpi_size() const noexcept;

  /// Returns the number of quadrants on all trees on this processor.
  int local_num_quadrants() const noexcept;

  /// Returns the number of quadrants on all trees on all processors
  std::ptrdiff_t global_num_quadrants() const noexcept;

  /// Returns a view of all trees.
  span<const tree<2>> trees() const noexcept;

  p4est_t* native() noexcept;
  const p4est_t* native() const noexcept;

private:
  struct destroyer {
    void operator()(p4est_t* p) const noexcept {
      if (p) {
        p4est_destroy(p);
      }
    }
  };
  std::unique_ptr<p4est_t, destroyer> m_handle{nullptr};
};

template <typename Predicate>
void refine_if(forest<2>& forest, Predicate predicate) {
  forest.native()->user_pointer = &predicate;
  p4est_refine(
      forest.native(), 0,
      [](p4est_t* forest, int which_tree, p4est_quadrant_t* quad) -> int {
        Predicate* pred = reinterpret_cast<Predicate*>(forest->user_pointer);
        quadrant<2> q(*quad);
        return fub::invoke(*pred, which_tree, q);
      },
      nullptr);
}

void balance(forest<2>& forest) noexcept;

optional<std::ptrdiff_t> find(const forest<2>& forest, int treeidx,
                              const quadrant<2>& quad) noexcept;

} // namespace p4est
} // namespace v1
} // namespace fub

#endif