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

#include "fub/p4est/tree.hpp"

namespace fub {
inline namespace v1 {
namespace p4est {
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

/// Returns a view of all locally stored quadrants.
span<const quadrant<2>> tree<2>::quadrants() const noexcept {
  return as_quadrant2_span(const_cast<sc_array*>(&m_native.quadrants));
}

/// Returns the highest local quadrant refinement level.
int tree<2>::max_level() const noexcept { return m_native.maxlevel; }

/// Returns the cumulative sum over earlier trees on this processor (locals
/// only).
std::ptrdiff_t tree<2>::offset() const noexcept {
  return m_native.quadrants_offset;
}

} // namespace p4est
} // namespace v1
} // namespace fub