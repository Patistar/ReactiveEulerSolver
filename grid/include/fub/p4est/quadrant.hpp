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

extern "C" {
#include <p4est.h>
#include <p4est_bits.h>
#include <p8est.h>
}

#include <array>

namespace fub {
inline namespace v1 {
namespace p4est {

/// \ingroup p4est
/// This is a wrapper for the quadrant types of p4est.
template <int Rank> class quadrant;

/// \ingroup p4est
/// This is a wrapper for the 2-dimensional quadrant type `p4est_quadrant_t`.
///
/// We use this type conveniently to index into out local patch data vector
/// of \a grid. This makes a quadrant being a handle type.
template <> class quadrant<2> {
public:
  /// \name Constructor

  quadrant() = default;

  /// Implicit conversion from its native type.
  constexpr quadrant(const p4est_quadrant_t& q) noexcept : m_handle{q} {}

  /// \name Observers

  /// Returns 2.
  static constexpr int rank() noexcept { return 2; }

  /// @{
  /// Returns the native type.
  ///
  /// \throws Nothing.
  constexpr const p4est_quadrant_t& get_native() const noexcept {
    return m_handle;
  }
  constexpr p4est_quadrant_t& get_native() noexcept { return m_handle; }
  /// @}

  /// Returns the index of the local patch data in the grid where this quadrant
  /// belongs to.
  constexpr long index() const noexcept { return m_handle.p.user_long; }

  /// Returns the tree index of this quadrant (only valid for ghost quadrants).
  constexpr int which_tree() const noexcept { return m_handle.p.which_tree; }

  /// Returns the owner rank of this quadrant (only valid for ghost quadrants).
  constexpr int local_num() const noexcept {
    return m_handle.p.piggy3.local_num;
  }

  /// Returns the refinement level of this quadrant.
  constexpr int get_level() const noexcept { return m_handle.level; }

  /// Returns a array of {x, y} coordinates which can be used to adapt
  /// coordinates.
  std::array<std::uint64_t, 2> get_coordinates() const noexcept {
    assert(m_handle.x >= 0 && m_handle.y >= 0);
    std::uint64_t x = reverse(reverse(m_handle.x, 30), get_level());
    std::uint64_t y = reverse(reverse(m_handle.y, 30), get_level());
    return {x, y};
  }

private:
  static uint32_t reverse(uint32_t x, int bits) noexcept {
    x = ((x & 0x55555555) << 1) | ((x & 0xAAAAAAAA) >> 1);
    x = ((x & 0x33333333) << 2) | ((x & 0xCCCCCCCC) >> 2);
    x = ((x & 0x0F0F0F0F) << 4) | ((x & 0xF0F0F0F0) >> 4);
    x = ((x & 0x00FF00FF) << 8) | ((x & 0xFF00FF00) >> 8);
    x = ((x & 0x0000FFFF) << 16) | ((x & 0xFFFF0000) >> 16);
    return x >> (32 - bits);
  }

  p4est_quadrant_t m_handle;
};

inline quadrant<2> face_neighbor(quadrant<2> quad, int face) noexcept {
  quadrant<2> nb;
  p4est_quadrant_face_neighbor(&quad.get_native(), face, &nb.get_native());
  return nb;
}

/// \ingroup p4est
/// This is a wrapper for the 3-dimensional quadrant type `p8est_quadrant_t`.
///
/// We use this type conveniently to index into out local patch data vector
/// of \a grid. This makes a quadrant being a handle type.
template <> class quadrant<3> {
public:
public:
  /// \name Constructor

  quadrant() = default;

  /// Implicit conversion from its native type.
  constexpr quadrant(const p8est_quadrant_t& q) noexcept : m_handle{q} {}

  /// \name Observers

  /// Returns 3.
  static constexpr int rank() noexcept { return 3; }

  /// Returns the native type.
  constexpr const p8est_quadrant_t& get_native() const noexcept {
    return m_handle;
  }

  /// Returns the index of the local patch data in the grid where this quadrant
  /// belongs to.
  constexpr long index() const noexcept { return m_handle.p.user_long; }

  /// Returns the refinement level of this quadrant.
  constexpr int get_level() const noexcept { return m_handle.level; }

  /// Returns a array of {x, y} coordinates which can be used to adapt
  /// coordinates.
  std::array<std::uint64_t, 2> get_coordinates() const noexcept {
    assert(m_handle.x >= 0 && m_handle.y >= 0);
    std::uint64_t x = reverse(reverse(m_handle.x, 30), get_level());
    std::uint64_t y = reverse(reverse(m_handle.y, 30), get_level());
    return {x, y};
  }

private:
  static uint32_t reverse(uint32_t x, int bits) noexcept {
    x = ((x & 0x55555555) << 1) | ((x & 0xAAAAAAAA) >> 1);
    x = ((x & 0x33333333) << 2) | ((x & 0xCCCCCCCC) >> 2);
    x = ((x & 0x0F0F0F0F) << 4) | ((x & 0xF0F0F0F0) >> 4);
    x = ((x & 0x00FF00FF) << 8) | ((x & 0xFF00FF00) >> 8);
    x = ((x & 0x0000FFFF) << 16) | ((x & 0xFFFF0000) >> 16);
    return x >> (32 - bits);
  }

  p8est_quadrant_t m_handle;
};

} // namespace p4est
} // namespace v1
} // namespace fub

#endif