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

#include "fub/p4est/forest.hpp"

namespace fub {
inline namespace v1 {
namespace p4est {

// Constructors

forest<2>::forest(const forest& other)
    : m_handle{p4est_copy(const_cast<p4est_t*>(other.native()), 0)} {}

forest<2>& forest<2>::operator=(const forest& other) {
  m_handle.reset(p4est_copy(const_cast<p4est_t*>(other.native()), 1));
  return *this;
}

forest<2>::forest(p4est_t* pointer) noexcept : m_handle{pointer} {}

forest<2>::forest(MPI_Comm communicator,
                  const ::fub::p4est::connectivity<2>& conn) noexcept
    : m_handle{p4est_new(communicator,
                         const_cast<p4est_connectivity_t*>(conn.native()), 0,
                         nullptr, nullptr)} {}

forest<2>::forest(MPI_Comm communicator, ::fub::p4est::connectivity<2>& conn,
                  int min_quads, int min_level, int fill_uniform) noexcept
    : m_handle{p4est_new_ext(
          communicator, const_cast<p4est_connectivity_t*>(conn.native()),
          min_quads, min_level, fill_uniform, 0, nullptr, nullptr)} {}

// Member Accessors

MPI_Comm forest<2>::mpi_communicator() const noexcept {
  return m_handle->mpicomm;
}

int forest<2>::mpi_rank() const noexcept { return m_handle->mpirank; }

int forest<2>::mpi_size() const noexcept { return m_handle->mpisize; }

int forest<2>::local_num_quadrants() const noexcept {
  return m_handle->local_num_quadrants;
}

std::ptrdiff_t forest<2>::global_num_quadrants() const noexcept {
  return m_handle->global_num_quadrants;
}

span<const tree<2>> forest<2>::trees() const noexcept {
  const tree<2>* pointer = reinterpret_cast<const tree<2>*>(
      p4est_tree_array_index(m_handle->trees, 0));
  std::ptrdiff_t size = m_handle->trees->elem_count;
  return {pointer, size};
}

p4est_t* forest<2>::native() noexcept { return m_handle.get(); }

const p4est_t* forest<2>::native() const noexcept { return m_handle.get(); }

optional<std::ptrdiff_t> find(const forest<2>& forest,
                              const quadrant<2>& quad) noexcept {
  int treeidx = quad.which_tree();
  const tree<2>& tree = forest.trees()[treeidx];
  span<const quadrant<2>> quads = tree.quadrants();
  auto pos = std::lower_bound(quads.begin(), quads.end(), quad);
  if (pos == quads.end() || quad != *pos) {
    return nullopt;
  }
  return pos - quads.begin();
}

void balance(forest<2>& forest) noexcept {
  p4est_balance(forest.native(), P4EST_CONNECT_FACE, nullptr);
}

} // namespace p4est
} // namespace v1
} // namespace fub