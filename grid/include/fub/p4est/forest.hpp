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

extern "C" {
#include <p4est.h>
#include <p4est_extended.h>
#include <p4est_ghost.h>
#include <p8est.h>
}

#include <memory>

namespace fub {
inline namespace v1 {
namespace p4est {

template <int Rank> class forest;

template <> class forest<2> {
public:
  forest() = default;
  explicit forest(p4est_t*);

  forest(MPI_Comm communicator, p4est_connectivity_t* conn) noexcept
      : m_handle{p4est_new(communicator, conn, 0, nullptr, nullptr)} {}

  template <typename Init>
  forest(MPI_Comm communicator, p4est_connectivity_t* conn, Init init) noexcept
      : m_handle{p4est_new(
            communicator, conn, 0,
            [](p4est_t* g, p4est_topidx_t which_tree, p4est_quadrant_t* quad) {
              return fub::invoke(*static_cast<Init*>(g->user_pointer), g,
                                 which_tree, quad);
            },
            &init)} {}

  template <typename Init>
  forest(MPI_Comm communicator, p4est_connectivity_t* conn,
         p4est_locidx_t min_quads, int min_level, int fill_uniform,
         Init init) noexcept
      : m_handle{p4est_new_ext(
            communicator, conn, min_quads, min_level, fill_uniform, 0,
            [](p4est_t* g, p4est_topidx_t which_tree, p4est_quadrant_t* quad) {
              return fub::invoke(*static_cast<Init*>(g->user_pointer), g,
                                 which_tree, quad);
            },
            &init)} {}

  p4est_t* operator->() const noexcept { return m_handle.get(); }

  operator p4est_t*() const noexcept { return m_handle.get(); }

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

inline p4est_locidx_t face_quadrant_exists(p4est_t* p4est, p4est_ghost_t* ghost,
                                           p4est_topidx_t treeid,
                                           const p4est_quadrant_t* q, int* face,
                                           int* hang, int* owner_rank) {
  return p4est_face_quadrant_exists(p4est, ghost, treeid, q, face, hang,
                                    owner_rank);
}

template <> class forest<3> {
public:
  forest() = default;
  explicit forest(p8est_t*);

  operator p8est_t*() noexcept;
  operator const p8est_t*() const noexcept;

private:
  struct destroyer {
    void operator()(p8est_t*) const noexcept;
  };
  std::unique_ptr<p8est_t, destroyer> m_handle{nullptr};
};

} // namespace p4est
} // namespace v1
} // namespace fub

#endif