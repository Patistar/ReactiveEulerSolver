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

#ifndef FUB_OCTREE_HPP
#define FUB_OCTREE_HPP

#include "fub/array.hpp"
#include "fub/face.hpp"
#include "fub/optional.hpp"

#include <algorithm>
#include <cassert>
#include <climits>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <type_traits>

namespace fub {

// Decleration                                          [class.basic_octant] {{{

/// @brief This type wraps octant coordinates in a morton index.
///
/// This class holds a valid (depth, morton_index) pair: Only those bits are
/// non-zero which describe indices for the specified octree depth.
///
/// We also transport the rank information of the problem through the template
/// parameter, which controls how many coordinates are interleaved in the morton
/// index.
template <typename Integral, int Rank> class basic_octant {
  static_assert(std::is_integral<Integral>::value,
                "Integral is not an integral type.");

public:
  static constexpr int total_bits = sizeof(Integral) * CHAR_BIT;
  static constexpr int max_depth = total_bits / Rank;

private:
  int m_depth;
  Integral m_morton_index;

public:
  // Constructors

  /// @brief Constructs a morton index and depth pointing to the root.
  constexpr basic_octant() = default;

  /// @brief Constructs a morton index pointing to the octact with specified
  /// depth and coordinates.
  constexpr basic_octant(int depth, const array<Integral, Rank> &coordinates);

  // Accessors

  /// @brief Returns the depth.
  constexpr int depth() const noexcept;

  /// @brief Returns its morton index.
  constexpr Integral morton_index() const noexcept;

  // Observers, these are more performant as member functions as they can modify
  // the morton index directly.

  /// @brief Returns the parent octant.
  constexpr basic_octant parent() const noexcept;

  /// @brief Returns the child octant at the specified child id.
  constexpr basic_octant child(int child_id) const;

  /// @brief Returns the lower descendant with specified depth.
  constexpr basic_octant lower_descendant(int depth) const;

  /// @brief Returns the upper descendant with specified depth.
  constexpr basic_octant upper_descendant(int depth) const;

  /// @brief Returns the next octant at this depth.
  ///
  /// The next octant is defined as the lowest octant `o` such that
  ///
  ///     *this < o        and         this->depth() == o.depth().
  ///
  /// If no such octant exists an exception is being thrown.
  constexpr basic_octant next() const;

  static constexpr basic_octant root() noexcept { return {}; }
};

template <int Rank> using octant = basic_octant<std::uint64_t, Rank>;

#define FUB_DEFINE_OCTANT_OPERATOR(op)                                         \
  template <typename Integral, int Rank>                                       \
  constexpr bool operator op(                                                  \
      const basic_octant<Integral, Rank> &o1,                                  \
      const basic_octant<Integral, Rank> &o2) noexcept {                       \
    return o1.morton_index() op o2.morton_index();                             \
  }

FUB_DEFINE_OCTANT_OPERATOR(<);
FUB_DEFINE_OCTANT_OPERATOR(<=);
FUB_DEFINE_OCTANT_OPERATOR(>);
FUB_DEFINE_OCTANT_OPERATOR(>=);

#undef FUB_DEFINE_OCTANT_OPERATOR

template <typename Integral, int Rank>
constexpr bool operator==(const basic_octant<Integral, Rank> &o1,
                          const basic_octant<Integral, Rank> &o2) noexcept {
  return o1.depth() == o2.depth() && o1.morton_index() == o2.morton_index();
}

template <typename Integral, int Rank>
constexpr bool operator!=(const basic_octant<Integral, Rank> &o1,
                          const basic_octant<Integral, Rank> &o2) noexcept {
  return !(o1 == o2);
}

// }}}

// Implementation                                       [class.basic_octant] {{{

// Helper Functions {{{
namespace octant_detail {

/// @brief Returns a bit mask containing as many 1's as Rank is.
///
/// An example:
///     octant_level_mask<1> = 0b1
///     octant_level_mask<2> = 0b11
///     octant_level_mask<3> = 0b111
template <typename Integral, int Rank>
struct level_mask
    : std::integral_constant<Integral, (Integral(1) << Rank) - 1> {};

template <typename Integral, std::size_t Rank>
Integral interleave_bits(int depth, Integral bit_position,
                         const array<Integral, Rank> &xs) noexcept {
  constexpr int max_depth = basic_octant<Integral, Rank>::max_depth;
  Integral interleaved{};
  for (std::size_t i = 0; i < Rank; ++i) {
    interleaved |= ((xs[i] & bit_position) >> depth)
                   << (Rank * (max_depth - 1) + i);
  }
  return interleaved;
}

} // namespace octant_detail
// }}}

// Constrcutor {{{
template <typename Integral, int Rank>
constexpr basic_octant<Integral, Rank>::basic_octant(
    int depth, const array<Integral, Rank> &coords)
    : m_depth{depth}, m_morton_index{} {
  // Check for valid input
  if (depth < 0 || max_depth < depth) {
    throw std::out_of_range("basic_octant: Depth is out of range.");
  }
  const Integral max = depth == 0 ? 0 : Integral(1) << max_depth - depth;
  if (!std::all_of(coords.begin(), coords.end(),
                   [max](Integral x) { return 0 <= x && x <= max; })) {
    throw std::out_of_range(
        "basic_octant: Coordinates are out of range for specified depth.");
  }
  // Construct the morton index by interleaving the coordinates bitwise.
  Integral bit_position(1);
  for (int d = 0; d < depth; ++d) {
    m_morton_index >>= Rank;
    m_morton_index |= octant_detail::interleave_bits(d, bit_position, coords);
    bit_position <<= 1;
  }
}
// }}}

// Accessors {{{
template <typename Integral, int Rank>
constexpr int basic_octant<Integral, Rank>::depth() const noexcept {
  return m_depth;
}

template <typename Integral, int Rank>
constexpr Integral basic_octant<Integral, Rank>::morton_index() const noexcept {
  return m_morton_index;
}
// }}}

// Observers {{{
template <int Dim, typename Integral, int Rank>
constexpr std::enable_if_t<(Dim < Rank), Integral>
coordinate(const basic_octant<Integral, Rank> &o) noexcept {
  if (o.depth() == 0) {
    return 0;
  }
  constexpr int max_depth = basic_octant<Integral, Rank>::max_depth;
  constexpr Integral first = Integral(1) << (Rank * (max_depth - 1) + Dim);
  const Integral morton_index = o.morton_index();
  Integral bit_position = first;
  Integral coordinate{};
  for (int depth = 0; depth < o.depth(); ++depth) {
    coordinate <<= 1;
    coordinate |= ((morton_index & bit_position) >>
                   (Rank * (max_depth - depth - 1) + Dim));
    bit_position >>= Rank;
  }
  return coordinate;
}

template <typename Integral, int Rank>
constexpr array<Integral, Rank>
coordinates(const basic_octant<Integral, Rank> &o) noexcept {
  constexpr int max_depth = basic_octant<Integral, Rank>::max_depth;
  constexpr Integral first = octant_detail::level_mask<Integral, Rank>::value
                             << (Rank * (max_depth - 1));
  const Integral morton_index = o.morton_index();
  Integral level_position = first;
  array<Integral, Rank> coordinates{};
  for (int depth = 0; depth < o.depth(); ++depth) {
    const Integral level =
        (morton_index & level_position) >> (Rank * (max_depth - depth - 1));
    for (int dim = 0; dim < Rank; ++dim) {
      coordinates[dim] <<= 1;
      coordinates[dim] |= ((level & (1 << dim)) >> dim);
    }
    level_position >>= Rank;
  }
  return coordinates;
}

template <typename Integral, int Rank>
constexpr int child_id(basic_octant<Integral, Rank> octant, int depth) {
  if (depth == 0) {
    return 0;
  } else if (!(0 < depth && depth < octant.max_depth)) {
    throw std::out_of_range("child_id: depth is out of range.");
  }
  constexpr Integral level = octant_detail::level_mask<Integral, Rank>::value;
  const int shift = Rank * (octant.max_depth - depth);
  assert(shift < octant.total_bits);
  const Integral mask = level << shift;
  const Integral masked = octant.morton_index() & mask;
  const Integral result = masked >> shift;
  return result;
}

template <typename Integral, int Rank>
constexpr basic_octant<Integral, Rank>
basic_octant<Integral, Rank>::parent() const noexcept {
  if (m_depth <= 1) {
    return basic_octant{};
  }
  const int depth = static_cast<int>(m_depth) - 1;
  const int shift = Rank * (max_depth - depth);
  assert(shift < total_bits);
  const Integral mask = ~((Integral(1) << shift) - 1);
  basic_octant o{};
  o.m_morton_index = m_morton_index & mask;
  o.m_depth = depth;
  return o;
}

template <typename Integral, int Rank>
constexpr basic_octant<Integral, Rank>
basic_octant<Integral, Rank>::child(int child_id) const {
  if (!(0 <= child_id &&
        child_id <= octant_detail::level_mask<int, Rank>::value)) {
    throw std::out_of_range("basic_octant::child: child_id is out of range.");
  }
  if (m_depth == max_depth) {
    return *this;
  }
  const int depth = m_depth + 1;
  const int shift = Rank * (max_depth - depth);
  const Integral mask = static_cast<Integral>(child_id) << shift;
  basic_octant o{};
  o.m_morton_index = m_morton_index | mask;
  o.m_depth = depth;
  return o;
}

template <typename Integral, int Rank>
constexpr basic_octant<Integral, Rank>
basic_octant<Integral, Rank>::lower_descendant(int depth) const {
  depth = std::max(m_depth, depth);
  if (depth > max_depth) {
    throw std::out_of_range(
        "basic_octant::lower_descendant: depth is out of range.");
  }
  basic_octant o{};
  o.m_morton_index = m_morton_index;
  o.m_depth = depth;
  return o;
}

template <typename Integral, int Rank>
constexpr basic_octant<Integral, Rank>
basic_octant<Integral, Rank>::upper_descendant(int depth) const {
  if (depth <= m_depth || m_depth == max_depth) {
    return *this;
  }
  if (depth > max_depth) {
    throw std::out_of_range(
        "basic_octant::upper_descendant: depth is out of range.");
  }
  basic_octant o{};
  constexpr int shift = total_bits % Rank;
  const Integral mask_upper = ~Integral{0} >> (Rank * m_depth + shift);
  const Integral mask_lower = ~(~Integral{0} >> (Rank * depth + shift));
  const Integral mask = mask_upper & mask_lower;
  o.m_morton_index = m_morton_index | mask;
  o.m_depth = depth;
  return o;
}

template <typename Integral, int Rank>
constexpr basic_octant<Integral, Rank>
basic_octant<Integral, Rank>::next() const {
  constexpr int shift = total_bits % Rank;
  const Integral last_index = ~(~Integral{0} >> (Rank * m_depth + shift));
  if (m_morton_index == last_index) {
    throw std::logic_error{"basic_octant::next: This octant has no successor."};
  }
  basic_octant o{*this};
  const Integral step = Integral(1) << (Rank * (max_depth - m_depth));
  o.m_morton_index += step;
  return o;
}

template <typename Integral, int Rank>
constexpr optional<basic_octant<Integral, Rank>>
face_neighbor(const basic_octant<Integral, Rank> &octant,
              const face &face) noexcept {
  array<Integral, Rank> xs = coordinates(octant);
  const int dim = as_int(face.dimension);
  const Integral max = (1 << octant.depth()) - 1;
  if ((xs[dim] == 0 && face.side == direction::left) ||
      (xs[dim] == max && face.side == direction::right)) {
    return {};
  }
  xs[dim] += sign(face.side);
  return basic_octant<Integral, Rank>(octant.depth(), xs);
}
// }}}

// }}}

// Decleration                                                [class.octree] {{{

/// @brief This class holds the local octree in a distributed process.
///
/// The local octree stores only leaf octants. It can perform local
/// transformations as refining and coarsening octants. It also supports a
/// logarithmic lookup and a way of iterating through all octants.
///
/// octants at the boundary need special treating for refine and coarsen
/// algorithms.
template <typename T, int Rank,
          typename Allocator = std::allocator<std::pair<const octant<Rank>, T>>>
class octree {
  int m_max_depth{octant<Rank>::max_depth};
  int m_depth{};
  std::map<octant<Rank>, T> m_octants;

public:
  using map_type = std::map<octant<Rank>, T>;
  using key_type = const octant<Rank>;
  using mapped_type = T;
  using value_type = std::pair<key_type, mapped_type>;
  using iterator = typename map_type::iterator;
  using const_iterator = typename map_type::const_iterator;

  /// @brief Constructs an empty octree.
  octree() = default;

  /// @brief Returns the number of local octants in the octree.
  std::size_t size() const noexcept;

  // Refine methods

  /// @brief Replaces the octant at the specified position `pos` by its children
  /// and call the transfer function to interpolate data from coarse to fine
  /// level.
  ///
  /// @return Returns the range of children which have been created in the
  /// process.
  template <typename TransferFunction>
  std::pair<const_iterator, const_iterator> refine(const_iterator pos,
                                                   TransferFunction transfer);

  // Coarsen methods

  /// @brief Replaces each family of octants with its parents in the range
  /// [first, last). Calls the transfer function to interpolate from fine to
  /// coarse level.
  ///
  /// @return Returns the last coarsened octant.
  template <typename TransferFunction>
  const_iterator coarsen(const_iterator first, const_iterator last,
                         TransferFunction transfer);

  // Iterators

  iterator begin() noexcept;
  const_iterator begin() const noexcept;
  const_iterator cbegin() const noexcept;

  iterator end() noexcept;
  const_iterator end() const noexcept;
  const_iterator cend() const noexcept;

  // Modifiers

  iterator insert(const_iterator hint, const value_type &value);
  iterator insert(const_iterator hint, value_type &&value);
  std::pair<iterator, bool> insert(const value_type &value);
  std::pair<iterator, bool> insert(value_type &&value);

  void erase(const_iterator pos);

  // Observers

  iterator find(const octant<Rank> &octant) noexcept;
  const_iterator find(const octant<Rank> &octant) const noexcept;

  iterator lower_bound(const octant<Rank> &octant) noexcept;
  const_iterator lower_bound(const octant<Rank> &octant) const noexcept;

  iterator upper_bound(const octant<Rank> &octant) noexcept;
  const_iterator upper_bound(const octant<Rank> &octant) const noexcept;
};

// }}}

// Implementation                                             [class.octree] {{{

template <typename T, int Rank, typename Allocator>
std::size_t octree<T, Rank, Allocator>::size() const noexcept {
  return m_octants.size();
}

template <typename T, int Rank, typename Allocator>
template <typename TransferFunction>
std::pair<typename octree<T, Rank, Allocator>::const_iterator,
          typename octree<T, Rank, Allocator>::const_iterator>
octree<T, Rank, Allocator>::refine(const_iterator pos,
                                   TransferFunction /* transfer */) {
  // dummy impl
  return {pos, pos};
}

template <typename T, int Rank, typename Allocator>
template <typename TransferFunction>
typename octree<T, Rank, Allocator>::const_iterator
octree<T, Rank, Allocator>::coarsen(const_iterator first,
                                    const_iterator /* last */,
                                    TransferFunction /* transfer */) {
  // dummy impl
  return first;
}

template <typename T, int Rank, typename Allocator>
typename octree<T, Rank, Allocator>::iterator
octree<T, Rank, Allocator>::begin() noexcept {
  return m_octants.begin();
}

template <typename T, int Rank, typename Allocator>
typename octree<T, Rank, Allocator>::const_iterator
octree<T, Rank, Allocator>::begin() const noexcept {
  return m_octants.begin();
}

template <typename T, int Rank, typename Allocator>
typename octree<T, Rank, Allocator>::const_iterator
octree<T, Rank, Allocator>::cbegin() const noexcept {
  return m_octants.cbegin();
}

template <typename T, int Rank, typename Allocator>
typename octree<T, Rank, Allocator>::iterator
octree<T, Rank, Allocator>::end() noexcept {
  return m_octants.end();
}

template <typename T, int Rank, typename Allocator>
typename octree<T, Rank, Allocator>::const_iterator
octree<T, Rank, Allocator>::end() const noexcept {
  return m_octants.end();
}

template <typename T, int Rank, typename Allocator>
typename octree<T, Rank, Allocator>::const_iterator
octree<T, Rank, Allocator>::cend() const noexcept {
  return m_octants.cend();
}

// Modifiers

template <typename T, int Rank, typename Allocator>
typename octree<T, Rank, Allocator>::iterator
octree<T, Rank, Allocator>::insert(const_iterator hint,
                                   const value_type &value) {
  return m_octants.insert(hint, value);
}

template <typename T, int Rank, typename Allocator>
typename octree<T, Rank, Allocator>::iterator
octree<T, Rank, Allocator>::insert(const_iterator hint, value_type &&value) {
  return m_octants.insert(hint, std::move(value));
}

template <typename T, int Rank, typename Allocator>
std::pair<typename octree<T, Rank, Allocator>::iterator, bool>
octree<T, Rank, Allocator>::insert(const value_type &value) {
  return m_octants.insert(value);
}

template <typename T, int Rank, typename Allocator>
std::pair<typename octree<T, Rank, Allocator>::iterator, bool>
octree<T, Rank, Allocator>::insert(value_type &&value) {
  return m_octants.insert(std::move(value));
}

template <typename T, int Rank, typename Allocator>
void octree<T, Rank, Allocator>::erase(const_iterator pos) {
  return m_octants.erase(pos);
}

// Observers

template <typename T, int Rank, typename Allocator>
typename octree<T, Rank, Allocator>::iterator
octree<T, Rank, Allocator>::find(const octant<Rank> &octant) noexcept {
  return m_octants.find(octant);
}

template <typename T, int Rank, typename Allocator>
typename octree<T, Rank, Allocator>::const_iterator
octree<T, Rank, Allocator>::find(const octant<Rank> &octant) const noexcept {
  return m_octants.find(octant);
}

template <typename T, int Rank, typename Allocator>
typename octree<T, Rank, Allocator>::iterator
octree<T, Rank, Allocator>::lower_bound(const octant<Rank> &octant) noexcept {
  return m_octants.lower_bound(octant);
}

template <typename T, int Rank, typename Allocator>
typename octree<T, Rank, Allocator>::const_iterator
octree<T, Rank, Allocator>::lower_bound(const octant<Rank> &octant) const
    noexcept {
  return m_octants.lower_bound(octant);
}

template <typename T, int Rank, typename Allocator>
typename octree<T, Rank, Allocator>::iterator
octree<T, Rank, Allocator>::upper_bound(const octant<Rank> &octant) noexcept {
  return m_octants.upper_bound(octant);
}

template <typename T, int Rank, typename Allocator>
typename octree<T, Rank, Allocator>::const_iterator
octree<T, Rank, Allocator>::upper_bound(const octant<Rank> &octant) const
    noexcept {
  return m_octants.upper_bound(octant);
}

// }}}

} // namespace fub

#endif // !FUB_OCTREE_HPP
