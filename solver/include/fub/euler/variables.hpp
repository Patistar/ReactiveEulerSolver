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

#ifndef FUB_EULER_VARIABLES_HPP
#define FUB_EULER_VARIABLES_HPP

#include "fub/variable.hpp"

namespace fub {
inline namespace v1 {
namespace euler {
struct Density : scalar_variable {
  static const char* name(int) noexcept { return "Density"; }
};

template <int Rank> struct Velocity : vector_variable<Rank> {
  static const char* name(int dim) noexcept {
    static const char* names[3]{"VelocityX", "VelocityY", "VelocityZ"};
    assert(0 <= dim && dim < Rank);
    return names[dim];
  }
};

template <int Rank> struct Momentum : vector_variable<Rank> {
  static const char* name(int dim) noexcept {
    static const char* names[3]{"MomentumX", "MomentumY", "MomentumZ"};
    assert(0 <= dim && dim < Rank);
    return names[dim];
  }
};

struct Pressure : scalar_variable {
  static const char* name(int) noexcept { return "Pressure"; }
};

struct Energy : scalar_variable {
  static const char* name(int) noexcept { return "EnergyStagnationDensity"; }
};

struct Temperature : scalar_variable {
  static const char* name(int) noexcept { return "Temperature"; }
};

struct SpeedOfSound : scalar_variable {
  static const char* name(int) noexcept { return "VelocitySound"; }
};

struct HeatCapacityAtConstantP : scalar_variable {
  static const char* name(int) noexcept { return "SpecificHeatPressure"; }
};

struct Enthalpy : scalar_variable {
  static const char* name(int) noexcept { return "Enthalpy"; }
};

struct Species : vector_variable<dynamic_extent> {
  using vector_variable<dynamic_extent>::vector_variable;
};

template <int Rank, int Dim = 0>
static constexpr const tag_t<Momentum<Rank>, Dim> momentum{};

template <int Rank> constexpr auto dynamic_momentum(int dim);

template <> constexpr auto dynamic_momentum<1>(int dim) {
  return variant<tag_t<Momentum<1>, 0>>{};
}

template <> constexpr auto dynamic_momentum<2>(int dim) {
  variant<tag_t<Momentum<2>, 0>, tag_t<Momentum<2>, 1>> v{};
  if (dim == 1) {
    v = momentum<2, 1>;
  }
  return v;
}

template <> constexpr auto dynamic_momentum<3>(int dim) {
  variant<tag_t<Momentum<3>, 0>, tag_t<Momentum<3>, 1>, tag_t<Momentum<3>, 2>> v{};
  if (dim == 1) {
    v = momentum<3, 1>;
  }
  if (dim == 2) {
    v = momentum<3, 2>;
  }
  return v;
}

template <int Rank, int Dim = 0>
static constexpr const tag_t<Velocity<Rank>, Dim> velocity{};

static constexpr const tag_t<Density> density{};
static constexpr const tag_t<Pressure> pressure{};
static constexpr const tag_t<Energy> energy{};
static constexpr const tag_t<Temperature> temperature{};
static constexpr const tag_t<SpeedOfSound> speed_of_sound{};
static constexpr const tag_t<HeatCapacityAtConstantP> cp{};
static constexpr const tag_t<Enthalpy> enthalpy{};

} // namespace euler
} // namespace v1
} // namespace fub

#endif
