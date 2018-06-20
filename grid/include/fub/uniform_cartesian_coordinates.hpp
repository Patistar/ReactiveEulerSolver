// Copyright (c) 2017 Maikel Nadolski
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

#ifndef FUB_UNIFORM_CARTESIAN_COORDINATES_HPP
#define FUB_UNIFORM_CARTESIAN_COORDINATES_HPP

#include "fub/array.hpp"
#include "fub/octree.hpp"
#include "fub/tuple.hpp"
#include "fub/type_traits.hpp"
#include "fub/utility.hpp"

#include <fmt/format.h>

namespace fub {

template <int Rank> class uniform_cartesian_coordinates {
  array<index, Rank> m_extents;
  array<double, Rank> m_lower;
  array<double, Rank> m_upper;
  array<double, Rank> m_dx;

public:
  uniform_cartesian_coordinates() = default;

  constexpr uniform_cartesian_coordinates(const array<double, Rank>& lower,
                                          const array<double, Rank>& upper,
                                          const array<index, Rank>& extents)
      : m_extents{extents}, m_lower{lower}, m_upper{upper},
        m_dx{compute_dx(m_lower, m_upper, m_extents)} {}

  const array<double, Rank>& lower() const noexcept { return m_lower; }
  const array<double, Rank>& upper() const noexcept { return m_upper; }
  const array<index, Rank>& extents() const noexcept { return m_extents; }
  const array<double, Rank>& dx() const noexcept { return m_dx; }

  template <typename... Is>
  constexpr std::enable_if_t<sizeof...(Is) == Rank, array<double, Rank>>
  operator()(Is... indices) const {
    array<double, Rank> coordinates{};
    array<index, Rank> is{{static_cast<index>(indices)...}};
    check_index_range(is);
    array<double, Rank> lambda{};
    std::transform(m_extents.begin(), m_extents.end(), is.begin(),
                   lambda.begin(), [](index e, index i) {
                     assert(0 <= i && i <= e);
                     return static_cast<double>(i) / e;
                   });
    for (int i = 0; i < Rank; ++i) {
      // This formula guarantees that m_lower and m_upper will be precisely
      // computed.
      coordinates[i] = (1 - lambda[i]) * m_lower[i] + lambda[i] * m_upper[i];
    }
    return coordinates;
  }

private:
  void check_index_range(const array<index, Rank>& is) const {
    for (int dim = 0; dim < Rank; ++dim) {
      if (is[dim] < 0 || m_extents[dim] < is[dim]) {
        std::string what =
            fmt::format("uniform_cartesian_coordinates: index {} in dimension "
                        "{} is out of range, extent is {}.",
                        is[dim], dim, m_extents[dim]);
        throw std::out_of_range(what);
      }
    }
  }

  static array<double, Rank> compute_dx(const array<double, Rank>& lower,
                                        const array<double, Rank>& upper,
                                        const array<index, Rank>& extents) {
    array<double, Rank> dx;
    std::transform(lower.begin(), lower.end(), upper.begin(), dx.begin(),
                   [](double lo, double up) { return up - lo; });
    std::transform(dx.begin(), dx.end(), extents.begin(), dx.begin(),
                   [](double dx, index e) { return dx / e; });
    return dx;
  }

  template <typename Archive>
  friend void serialize(Archive& archive, uniform_cartesian_coordinates& coords,
                        unsigned) {
    archive & coords.m_extents;
    archive & coords.m_lower;
    archive & coords.m_upper;
    archive & coords.m_dx;
  }
};

template <int Rank>
uniform_cartesian_coordinates<Rank>
adapt(const uniform_cartesian_coordinates<Rank>& coords,
      const octant<Rank>& octant) {
  array<index, Rank> extents = coords.extents();
  index refinement_ratio = (1 << octant.depth());
  std::transform(extents.begin(), extents.end(), extents.begin(),
                 [=](index e) { return e * refinement_ratio; });
  uniform_cartesian_coordinates<Rank> intermediate(coords.lower(),
                                                   coords.upper(), extents);
  array<std::uint64_t, Rank> idx = coordinates(octant);
  std::transform(idx.begin(), idx.end(), coords.extents().begin(), idx.begin(),
                 [](index lo, index e) { return lo * e; });
  array<double, Rank> lower = fub::apply(intermediate, idx);
  std::transform(idx.begin(), idx.end(), coords.extents().begin(), idx.begin(),
                 [](index lo, index e) { return lo + e; });
  array<double, Rank> upper = fub::apply(intermediate, idx);
  return {lower, upper, coords.extents()};
}

} // namespace fub

#endif // !UNIFORMCARTESIANMAPPING_HPP
