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

#include "fub/uniform_cartesian_coordinates.hpp"

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

TEST_CASE("1D mapping is in range") {
  double lower = 0.2;
  double upper = 0.3;
  fub::uniform_cartesian_coordinates<1> coordinates({lower}, {upper}, {4});
  REQUIRE(coordinates(0)[0] == lower);
  double x = coordinates(1)[0];
  REQUIRE((lower < x && x < upper));
  x = coordinates(2)[0];
  REQUIRE((lower < x && x < upper));
  REQUIRE(coordinates(4)[0] == upper);
  REQUIRE_THROWS(coordinates(-1));
  REQUIRE_THROWS(coordinates(5));
}

TEST_CASE("2D ") {
  fub::array<double, 2> lower{};
  fub::array<double, 2> upper{{1.33, 2.51}};
  fub::uniform_cartesian_coordinates<2> coordinates(lower, upper, {4, 8});

  SECTION("mapping is in range") {
    REQUIRE(coordinates(0, 0) == lower);
    REQUIRE(coordinates(4, 8) == upper);
    REQUIRE_THROWS(coordinates(-1, 0));
    REQUIRE_THROWS(coordinates(5, 0));
  }

  SECTION("Refinement 2 with upper right") {
    fub::octant<2> octant(1, {1, 1});
    auto adapted = adapt(coordinates, octant);
    REQUIRE(adapted.upper() == coordinates.upper());
    REQUIRE(adapted.lower() == coordinates(2, 4));
    REQUIRE(adapted.extents() == coordinates.extents());
  }

  SECTION("Refinement 2 with lower left") {
    fub::octant<2> octant(1, {0, 0});
    auto adapted = adapt(coordinates, octant);
    REQUIRE(adapted.upper() == coordinates(2, 4));
    REQUIRE(adapted.lower() == coordinates.lower());
    REQUIRE(adapted.extents() == coordinates.extents());
  }

  SECTION("Refinement 2 with lower right") {
    fub::octant<2> octant(1, {1, 0});
    auto adapted = adapt(coordinates, octant);
    REQUIRE(adapted.upper() == coordinates(4, 4));
    REQUIRE(adapted.lower() == coordinates(2, 0));
    REQUIRE(adapted.extents() == coordinates.extents());
  }

  SECTION("Refinement 2 with upper left") {
    fub::octant<2> octant(1, {0, 1});
    auto adapted = adapt(coordinates, octant);
    REQUIRE(adapted.upper() == coordinates(2, 8));
    REQUIRE(adapted.lower() == coordinates(0, 4));
    REQUIRE(adapted.extents() == coordinates.extents());
  }
}

