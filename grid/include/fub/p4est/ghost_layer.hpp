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

#ifndef FUB_P4EST_GHOST_LAYER_HPP
#define FUB_P4EST_GHOST_LAYER_HPP

#include "fub/span.hpp"

#include "fub/p4est/quadrant.hpp"
#include "fub/p4est/forest.hpp"

#include <p4est.h>
#include <p4est_ghost.h>

#include <memory>

namespace fub {
inline namespace v1 {
namespace p4est {

template <int Rank> struct ghost_layer;

/// \ingroup p4est
/// This class is a wrapper for `p4est_ghost_t`.
template <> struct ghost_layer<2> {
public:
  /// Constructs a ghost_layer by claiming ownership of `pointer`.
  explicit ghost_layer(p4est_ghost_t* pointer);

  /// Constrcuts a new ghost layer for a specified forest.
  ///
  /// This constructor is a MPI synchronization point and must be passed by all
  /// MPI ranks.
  explicit ghost_layer(const forest<2>& forest);

  /// Returns the number of trees in the ghost layer
  int num_trees() const noexcept;

  /// Returns an array of quadrants which make up the ghost layer around p4est.
  span<const quadrant<2>> quadrants() const noexcept;

  /// Return an array of local quadrants that touch the parallel boundary from
  /// the inside, id est, that are ghosts in the perspective of at least one
  /// other processor.
  span<const quadrant<2>> mirrors() const noexcept;

  /// Returns a partition of `quadrants()` by tree equivalence.
  span<const int> quadrants_tree_offsets() const noexcept;

  /// Returns a partition of `quadrants()` by owner equivalence.
  span<const int> quadrants_process_offsets() const noexcept;

  /// Returns an array of indices for `mirrors()` which is grouped by mpi ranks.
  span<const int> mirrors_by_process() const noexcept;

  /// Returns a parititon for `mirrors_by_process()` with owner equivalence.
  span<const int> mirrors_by_process_offsets() const noexcept;

  p4est_ghost_t* native() noexcept { return m_handle.get(); }
  const p4est_ghost_t* native() const noexcept { return m_handle.get(); }

private:
  struct destroyer {
    void operator()(p4est_ghost_t* p) const noexcept;
  };
  std::unique_ptr<p4est_ghost_t, destroyer> m_handle{nullptr};
};

} // namespace p4est
} // namespace v1
} // namespace fub

#endif