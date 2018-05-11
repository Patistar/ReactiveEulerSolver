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

#include "fub/variables.hpp"

namespace fub {
namespace euler {
struct Density {
  static constexpr const char* name() noexcept { return "Density"; }
  static const std::array<double, 5> physical_dimensions() noexcept {
    // [kg / m^3]
    return {{1, -3}};
  }
};

template <int Dim = 0> struct Velocity;

template <> struct Velocity<0> {
  static constexpr const char* name() noexcept { return "VelocityX"; }
  static const std::array<double, 5> physical_dimensions() {
    // [m / s]
    return {{0, +1, -1}};
  }
};

template <> struct Velocity<1> {
  static constexpr const char* name() noexcept { return "VelocityY"; }
  static constexpr std::array<double, 5> physical_dimensions() {
    // [m / s]
    return {{0, +1, -1}};
  }
};

template <> struct Velocity<2> {
  static constexpr const char* name() noexcept { return "VelocityZ"; }
  static constexpr std::array<double, 5> physical_dimensions() {
    // [m / s]
    return {{0, +1, -1}};
  }
};

template <int Dim> struct Momentum;

template <> struct Momentum<0> {
  static constexpr const char* name() noexcept { return "MomentumX"; }
  static constexpr std::array<double, 5> physical_dimensions() {
    // [kg / m^2 s]
    return {{1, -2, -1}};
  }
};

template <> struct Momentum<1> {
  static constexpr const char* name() noexcept { return "MomentumY"; }
  static constexpr std::array<double, 5> physical_dimensions() {
    // [kg / m^2 s]
    return {{1, -2, -1}};
  }
};

template <> struct Momentum<2> {
  static constexpr const char* name() noexcept { return "MomentumZ"; }
  static constexpr std::array<double, 5> physical_dimensions() {
    // [kg / m^2 s]
    return {{1, -2, -1}};
  }
};

struct Pressure {
  static constexpr const char* name() noexcept { return "Pressure"; }
  static constexpr std::array<double, 5> physical_dimensions() {
    // [kg / m s^2]
    return {{1, -1, -2}};
  }
};

struct Energy {
  static constexpr const char* name() noexcept {
    return "EnergyStagnationDensity";
  }
};

struct Temperature {
  static constexpr const char* name() noexcept { return "Temperature"; }
};

struct SpeedOfSound {
  static constexpr const char* name() noexcept { return "VelocitySound"; }
};

struct HeatCapacityAtConstantP {
  static constexpr const char* name() noexcept {
    return "SpecificHeatPressure";
  }
};

struct Enthalpy {
  static constexpr const char* name() noexcept { return "Enthalpy"; }
};

namespace variables {
static constexpr const Density density{};
static constexpr const Pressure pressure{};
static constexpr const Energy energy{};
static constexpr const Temperature temperature{};
static constexpr const SpeedOfSound speed_of_sound{};
static constexpr const HeatCapacityAtConstantP cp{};
static constexpr const Enthalpy enthalpy{};

template <int Dim> static constexpr const Velocity<Dim> velocity{};
template <int Dim> static constexpr const Momentum<Dim> momentum{};

} // namespace variables
} // namespace euler
} // namespace fub

#endif
