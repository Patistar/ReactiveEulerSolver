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

#ifndef FUB_P4EST_QUADRANT_HPP
#define FUB_P4EST_QUADRANT_HPP

#include "fub/face.hpp"
#include "fub/optional.hpp"
#include "fub/span.hpp"

extern "C" {
#include <p4est.h>
#include <p4est_bits.h>
#include <p8est.h>
}

#include <array>

namespace fub {
inline namespace v1 {
namespace p4est {

template <int Rank> class quadrant;

/// \ingroup p4est
/// This is a wrapper for the 2-dimensional quadrant type `p4est_quadrant_t`.
template <> class quadrant<2> {
public:
  /// \name Constructors

  quadrant() = default;
  quadrant(const quadrant&) = default;

  quadrant(const p4est_quadrant_t& q) : m_native{q} {}

  /// \name Observers

  /// Returns 2.
  static constexpr int rank() noexcept { return 2; }

  /// Returns the tree index of this quadrant (only valid for ghost and mirror
  /// quadrants).
  int which_tree() const noexcept;

  /// Returns the local index with respect to some p4est tree of this quadrant.
  int local_num() const noexcept;

  /// Returns the refinement level of this quadrant.
  int level() const noexcept;

  /// Returns an `std::array<int, 2>` of x and y coordinates which can be used
  /// to adapt coordinates.
  std::array<int, 2> coordinates() const noexcept;

  /// \name Base Access

  /// Returns a reference to the internal const p4est_quadrant_t object.
  const p4est_quadrant_t& native() const noexcept { return m_native; }

private:
  p4est_quadrant_t m_native;
};

bool operator==(const quadrant<2>& lhs, const quadrant<2>& rhs) noexcept;
bool operator!=(const quadrant<2>& lhs, const quadrant<2>& rhs) noexcept;
bool operator<(const quadrant<2>& lhs, const quadrant<2>& rhs) noexcept;

quadrant<2> face_neighbor(const quadrant<2>& quad, int face) noexcept;

optional<face> find_adjacent_face(const quadrant<2>& left,
                                  const quadrant<2>& right) noexcept;

std::array<quadrant<2>, 4> children(const quadrant<2>& quad) noexcept;

template <int Rank>
span<const quadrant<Rank>> as_quadrant_span(sc_array_t& array) {
  static_assert(std::is_standard_layout<quadrant<Rank>>::value,
                "Can not cast p4est_quadrant_t* to quadrant<2>*");
  std::ptrdiff_t size = array.elem_count;
  const quadrant<Rank>* pointer =
      size == 0 ? nullptr
                : reinterpret_cast<const quadrant<Rank>*>(
                      p4est_quadrant_array_index(&array, 0));
  return {pointer, size};
}

} // namespace p4est
} // namespace v1
} // namespace fub

#endif