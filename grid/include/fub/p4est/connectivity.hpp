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

#ifndef FUB_P4EST_CONNECTIVITY_HPP
#define FUB_P4EST_CONNECTIVITY_HPP

#include <p4est.h>

#include <memory>

namespace fub {
inline namespace v1 {
namespace p4est {

struct unit_square_t {};
static constexpr unit_square_t unit_square{};

struct periodic_t {};
static constexpr periodic_t periodic{};

template <int Rank> class connectivity;

/// \ingroup p4est
/// This is a wrapper for `p4est_connectivity_t` which describes how multiple
/// "blocks" are connected to in 2D.
template <> class connectivity<2> {
public:
  /// \name Constructors & Assignment

  /// Constructs an invalid connectivity.
  connectivity() = default;

  /// The move constructor releases the ownership of the internal pointer.
  connectivity(connectivity&& other) = default;

  /// The move assignment releases the ownership of the internal pointer.
  connectivity& operator=(connectivity&&) = default;

  /// Constructs a single block with boundaries.
  ///
  /// \b Example:
  /// \code
  /// // This line constructs a connectivity structure by calling
  /// // p4est_connectivity_new_unitsquare()
  /// connectivity<2> conn(unit_square);
  /// \endcode
  ///
  /// \param[in] unit_square This tag is used to distinguish between the
  /// constructors.
  ///
  /// \throws Nothing.
  connectivity(unit_square_t) noexcept // NOLINT
      : m_handle{p4est_connectivity_new_unitsquare()} {}

  /// Constructs a single block with periodic boundaries.
  ///
  /// \b Example:
  /// \code
  /// // This line constructs a connectivity structure by calling
  /// // p4est_connectivity_new_periodic()
  /// connectivity<2> conn(periodic);
  /// \endcode
  ///
  /// \param[in] periodic This tag is used to distinguish between the
  /// constructors.
  ///
  /// \throws Nothing.
  connectivity(periodic_t) noexcept // NOLINT
      : m_handle{p4est_connectivity_new_periodic()} {}

  /// Returns a pointer to `const p4est_connectivity_t`.
  const p4est_connectivity_t* native_handle() const noexcept;

  /// Returns a pointer to `p4est_connectivity_t`.
  p4est_connectivity_t* native_handle() noexcept;

  /// Returns the number of trees.
  int num_trees() const noexcept;

private:
  struct destroyer {
    void operator()(p4est_connectivity_t* handle) const noexcept {
      p4est_connectivity_destroy(handle);
    }
  };

  /// This unique pointer manages the pointer and destroys the
  /// p4est_connectivity structure for us.
  std::unique_ptr<p4est_connectivity_t, destroyer> m_handle{nullptr};
};

} // namespace p4est
} // namespace v1
} // namespace fub

#endif