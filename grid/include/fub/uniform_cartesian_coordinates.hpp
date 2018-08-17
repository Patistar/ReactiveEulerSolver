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

#include "fub/tuple.hpp"
#include "fub/type_traits.hpp"
#include "fub/utility.hpp"

#include <array>
#include <fmt/format.h>

namespace fub {
inline namespace v1 {

template <int Rank> class uniform_cartesian_coordinates {
  std::array<index, Rank> m_extents;
  std::array<double, Rank> m_lower;
  std::array<double, Rank> m_upper;
  std::array<double, Rank> m_dx;

public:
  uniform_cartesian_coordinates() = default;

  constexpr uniform_cartesian_coordinates(
      const std::array<double, Rank>& lower,
      const std::array<double, Rank>& upper,
      const std::array<index, Rank>& extents)
      : m_extents{extents}, m_lower{lower}, m_upper{upper},
        m_dx{compute_dx(m_lower, m_upper, m_extents)} {}

  static constexpr int rank() noexcept { return Rank; }

  const std::array<double, Rank>& lower() const noexcept { return m_lower; }
  const std::array<double, Rank>& upper() const noexcept { return m_upper; }
  const std::array<index, Rank>& extents() const noexcept { return m_extents; }
  const std::array<double, Rank>& dx() const noexcept { return m_dx; }
  double dx(int dim) const noexcept { return m_dx[dim]; }

  template <typename... Is>
  constexpr std::enable_if_t<sizeof...(Is) == Rank, std::array<double, Rank>>
  operator()(Is... indices) const {
    std::array<double, Rank> coordinates{};
    std::array<index, Rank> is{{static_cast<index>(indices)...}};
    check_index_range(is);
    std::array<double, Rank> lambda{};
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

  template <typename... IndexTypes>
  constexpr std::array<double, Rank> cell(IndexTypes... is) const {
    std::array<std::ptrdiff_t, Rank> indices{is...};
    std::array<std::ptrdiff_t, Rank> next;
    std::transform(indices.begin(), indices.end(), next.begin(),
                   [](std::ptrdiff_t i) { return i + 1; });
    auto x0 = fub::apply(*this, indices);
    auto x1 = fub::apply(*this, next);
    std::array<double, Rank> x;
    std::transform(x0.begin(), x0.end(), x1.begin(), x.begin(),
                   [](double lhs, double rhs) { return 0.5 * (lhs + rhs); });
    return x;
  }

private:
  void check_index_range(const std::array<index, Rank>& is) const {
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

  static std::array<double, Rank>
  compute_dx(const std::array<double, Rank>& lower,
             const std::array<double, Rank>& upper,
             const std::array<index, Rank>& extents) {
    std::array<double, Rank> dx;
    std::transform(lower.begin(), lower.end(), upper.begin(), dx.begin(),
                   [](double lo, double up) { return up - lo; });
    std::transform(dx.begin(), dx.end(), extents.begin(), dx.begin(),
                   [](double dx, index e) { return dx / e; });
    return dx;
  }

  template <typename Archive>
  friend void serialize(Archive& archive, uniform_cartesian_coordinates& coords,
                        unsigned) {
    // clang-format off
    archive & coords.m_extents;
    archive & coords.m_lower;
    archive & coords.m_upper;
    archive & coords.m_dx;
    // clang-format on
  }
};

template <int Rank, typename Oct>
uniform_cartesian_coordinates<Rank>
adapt(const uniform_cartesian_coordinates<Rank>& coords,
      const Oct& octant) {
  std::array<index, Rank> extents = coords.extents();
  index refinement_ratio = (1 << octant.level());
  std::transform(extents.begin(), extents.end(), extents.begin(),
                 [=](index e) { return e * refinement_ratio; });
  uniform_cartesian_coordinates<Rank> intermediate(coords.lower(),
                                                   coords.upper(), extents);
  std::array<int, Rank> idx = octant.coordinates();
  std::transform(idx.begin(), idx.end(), coords.extents().begin(), idx.begin(),
                 [](index lo, index e) { return lo * e; });
  std::array<double, Rank> lower = fub::apply(intermediate, idx);
  std::transform(idx.begin(), idx.end(), coords.extents().begin(), idx.begin(),
                 [](index lo, index e) { return lo + e; });
  std::array<double, Rank> upper = fub::apply(intermediate, idx);
  return {lower, upper, coords.extents()};
}

} // namespace v1
} // namespace fub

#endif // !UNIFORMCARTESIANMAPPING_HPP
