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

#include "fub/p4est/ghost_layer.hpp"

namespace fub {
inline namespace v1 {
namespace p4est {

ghost_layer<2>::ghost_layer(p4est_ghost_t* pointer) : m_handle{pointer} {}

ghost_layer<2>::ghost_layer(const forest<2>& forest)
    : m_handle{p4est_ghost_new(const_cast<p4est_t*>(forest.native()),
                               P4EST_CONNECT_FACE)} {}

int ghost_layer<2>::num_trees() const noexcept { return m_handle->num_trees; }

namespace {
span<const quadrant<2>> as_quadrant2_span(sc_array* array) noexcept {
  std::ptrdiff_t size = array->elem_count;
  if (size) {
    const quadrant<2>* pointer = reinterpret_cast<const quadrant<2>*>(
        p4est_quadrant_array_index(array, 0));
    return {pointer, size};
  }
  return {};
}
} // namespace

span<const quadrant<2>> ghost_layer<2>::quadrants() const noexcept {
  return as_quadrant2_span(&m_handle->ghosts);
}

span<const quadrant<2>> ghost_layer<2>::mirrors() const noexcept {
  return as_quadrant2_span(&m_handle->mirrors);
}

span<const int> ghost_layer<2>::quadrants_tree_offsets() const noexcept {
  return {m_handle->tree_offsets, num_trees() + 1};
}

span<const int> ghost_layer<2>::quadrants_process_offsets() const noexcept {
  return {m_handle->proc_offsets, m_handle->mpisize + 1};
}

span<const int> ghost_layer<2>::mirrors_by_process_offsets() const noexcept {
  return {m_handle->mirror_proc_offsets, m_handle->mpisize + 1};
}

span<const int> ghost_layer<2>::mirrors_by_process() const noexcept {
  return {m_handle->mirror_proc_mirrors, mirrors_by_process_offsets()[m_handle->mpisize]};
}

void ghost_layer<2>::destroyer::operator()(p4est_ghost_t* p) const noexcept {
  p4est_ghost_destroy(p);
}

} // namespace p4est
} // namespace v1
} // namespace fub
