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

#ifndef FUB_P4EST_TREE_HPP
#define FUB_P4EST_TREE_HPP

#include "fub/p4est/quadrant.hpp"
#include "fub/span.hpp"

namespace fub {
inline namespace v1 {
namespace p4est {

template <int Rank> class tree;

/// \ingroup p4est
/// This is a wrapper for the 2-dimensional tree type `p4est_tree_t`.
template <> class tree<2> {
public:
  /// \name Observers

  tree() = default;
  tree(const tree&) = default;

  /// The explicit conversion operator from the native type.
  explicit tree(const p4est_tree_t& q) : m_native{q} {}

  /// Returns 2.
  static constexpr int rank() noexcept { return 2; }

  /// Returns a view of all locally stored quadrants.
  span<const quadrant<2>> quadrants() const noexcept;

  /// Returns the highest local quadrant refinement level.
  int max_level() const noexcept;

  /// Returns the cumulative sum over earlier trees on this processor (locals only).
  std::ptrdiff_t offset() const noexcept;


private:
  p4est_tree_t m_native;
};

} // namespace p4est
} // namespace v1
} // namespace fub

#endif