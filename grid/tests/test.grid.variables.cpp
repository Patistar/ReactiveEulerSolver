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

#include "fub/variables.hpp"

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

using fub::quantities;
using fub::vector_variable;

struct Density {
  using value_type = double;
};
static constexpr Density density;

struct Velocity {
  using value_type = double;
};
static constexpr Velocity velocity;

struct Vector : vector_variable<Density, Velocity> {};

TEST_CASE("Can be constructed") {
  quantities<> empty{};
  quantities<Density> v1{42};
  REQUIRE(v1[density] == 42);
  quantities<Velocity> v2{24};
  REQUIRE(v2[velocity] == 24);
  quantities<Density, Velocity> rhou1;
  rhou1[density] = 0;
  rhou1[velocity] = 1;
  REQUIRE(rhou1[density] == 0);
  REQUIRE(rhou1[velocity] == 1);
  quantities<Density, Velocity> rhou2{};
  REQUIRE(rhou2[density] == 0);
  REQUIRE(rhou2[velocity] == 0);
  quantities<Velocity, Density> urho{42, 24};
  REQUIRE(urho[density] == 24);
  REQUIRE(urho[velocity] == 42);

  quantities<Density, Velocity> rhou3 = urho;
  REQUIRE(rhou3[density] == urho[density]);
  REQUIRE(rhou3[velocity] == urho[velocity]);

  quantities<Vector> rhou4 = urho;
  REQUIRE(rhou4[density] == urho[density]);
  REQUIRE(rhou4[velocity] == urho[velocity]);
}

TEST_CASE("Accessing flux vector variables") {
  fub::add_flux_t<quantities<Vector>> f;
  f[density] = 0;
  f[velocity] = 0;
}

TEST_CASE("Accessing simd flux vector variables") {
  fub::add_flux_t<fub::add_simd_t<quantities<Vector>>> f;
  f[density] = 0;
  f[velocity] = 0;
}
