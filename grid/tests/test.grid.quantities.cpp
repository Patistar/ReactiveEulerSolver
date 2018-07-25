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

#include "fub/quantities.hpp"

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

using fub::quantities;
using fub::types;

struct Density {
  using value_type = double;
};
static constexpr Density density;

struct Velocity {
  using value_type = double;
};
static constexpr Velocity velocity;

struct Vector : fub::variable_union<Density, Velocity> {};

TEST_CASE("Can be constructed") {
  quantities<Density> v1{42};
  REQUIRE(v1[density] == 42);
  quantities<Velocity> v2{24};
  REQUIRE(v2[velocity] == 24);
  quantities<Density, Velocity> rhou1;
  rhou1[density] = 0;
  rhou1[velocity] = 1;
  REQUIRE(rhou1[density] == 0);
  REQUIRE(rhou1[velocity] == 1);
}

TEST_CASE("Accessing simd vector variables") {
  fub::simd_quantities<types<Vector>> f;
  fub::simd<double> check_rho = 42.0;
  fub::simd<double> check_u = 24.0;
  SECTION("assign simd") {
    f[density] = check_rho;
    f[velocity] = check_u;
    REQUIRE(all_of(f[density] == check_rho));
    REQUIRE(all_of(f[velocity] == check_u));
  }
  SECTION("broadcast double") {
    f[density] = 42.0;
    f[velocity] = 24.0;
    REQUIRE(all_of(f[density] == check_rho));
    REQUIRE(all_of(f[velocity] == check_u));
  }
}


TEST_CASE("Accessing simd vector variables 2") {
  fub::simd_quantities<types<Density, Velocity>> f{};
  fub::simd<double> check_rho = 42.0;
  fub::simd<double> check_u = 24.0;
  REQUIRE(all_of(f[density] == fub::simd<double>{}));
  REQUIRE(all_of(f[velocity] == fub::simd<double>{}));
  SECTION("assign simd") {
    f[density] = check_rho;
    f[velocity] = check_u;
    REQUIRE(all_of(f[density] == check_rho));
    REQUIRE(all_of(f[velocity] == check_u));
  }
  SECTION("broadcast double") {
    f[density] = 42.0;
    f[velocity] = 24.0;
    REQUIRE(all_of(f[density] == check_rho));
    REQUIRE(all_of(f[velocity] == check_u));
  }
}

// TEST_CASE("Accessing simd flux vector variables") {
//   fub::add_flux_t<fub::add_simd_t<quantities<Vector>>> f;
//   f[density] = 0;
//   f[velocity] = 0;
// }
