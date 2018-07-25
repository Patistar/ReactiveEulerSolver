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

#include "fub/quantities_ref.hpp"

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

using fub::quantities;
using fub::quantities_ref;
using fub::types;

struct Density {
  using value_type = double;
};
static constexpr Density density;

struct Velocity {
  using value_type = double;
};
static constexpr Velocity velocity;

struct Pressure {
  using value_type = double;
};
static constexpr Pressure pressure;

struct Vector : fub::variable_union<Density, Velocity> {};

TEST_CASE("Can be constructed from mutable") {
  quantities<Density> q(42);
  quantities_ref<Density> ref(q);
  REQUIRE(&ref[density] == &q[density]);
  REQUIRE(sizeof(ref) == sizeof(void*));
  REQUIRE(q[density] == 42);
  REQUIRE(ref[density] == 42);
  ref[density] = 24;
  REQUIRE(q[density] == 24);
  REQUIRE(ref[density] == 24);
}

TEST_CASE("Can be constructed from const") {
  const quantities<Density> q(42);
  quantities_ref<const Density> ref(q);
  REQUIRE(&ref[density] == &q[density]);
  REQUIRE(sizeof(ref) == sizeof(void*));
  REQUIRE(q[density] == 42);
  REQUIRE(ref[density] == 42);
}

TEST_CASE("Can be constructed from mutable but to const") {
  quantities<Density> q(42);
  quantities_ref<const Density> ref(q);
  REQUIRE(&ref[density] == &q[density]);
  REQUIRE(sizeof(ref) == sizeof(void*));
  REQUIRE(q[density] == 42);
  REQUIRE(ref[density] == 42);
  q[density] = 24;
  REQUIRE(q[density] == 24);
  REQUIRE(ref[density] == 24);
}

TEST_CASE("Can be constructed from mutable Vector") {
  quantities<Vector> q({{42., 421.}});
  quantities_ref<Vector> ref(q);
  REQUIRE(&ref[density] == &q[density]);
  REQUIRE(sizeof(ref) == sizeof(void*));
  REQUIRE(q[density] == 42);
  REQUIRE(ref[density] == 42);
  ref[density] = 24;
  REQUIRE(q[density] == 24);
  REQUIRE(ref[density] == 24);
}

TEST_CASE("Can be constructed from const Vector") {
  const quantities<Vector> q({{42., 421.}});
  quantities_ref<const Vector> ref(q);
  REQUIRE(&ref[density] == &q[density]);
  REQUIRE(sizeof(ref) == sizeof(void*));
  REQUIRE(q[density] == 42);
  REQUIRE(ref[density] == 42);
}

TEST_CASE("Can be constructed from mutable but to const Vector") {
  quantities<Vector> q({{42., 421.}});
  quantities_ref<const Vector> ref(q);
  REQUIRE(&ref[density] == &q[density]);
  REQUIRE(sizeof(ref) == sizeof(void*));
  REQUIRE(q[density] == 42);
  REQUIRE(ref[density] == 42);
  q[density] = 24;
  REQUIRE(q[density] == 24);
  REQUIRE(ref[density] == 24);
}

TEST_CASE("Can be constructed from mutable Vector, Pressure") {
  quantities<Vector, Pressure> q{};
  quantities_ref<Vector, Pressure> ref(q);
  REQUIRE(&ref[density] == &q[density]);
  REQUIRE(sizeof(ref) == 2 * sizeof(void*));
  REQUIRE(q[density] == 0);
  REQUIRE(ref[density] == 0);
  ref[density] = 24;
  REQUIRE(q[density] == 24);
  REQUIRE(ref[density] == 24);
}

TEST_CASE("Can be constructed from const Vector, Pressure") {
  const quantities<Vector, Pressure> q{};
  quantities_ref<const Vector, const Pressure> ref(q);
  REQUIRE(&ref[density] == &q[density]);
  REQUIRE(sizeof(ref) == 2 * sizeof(void*));
  REQUIRE(q[density] == 0);
  REQUIRE(ref[density] == 0);
}

TEST_CASE("Can be constructed from mutable but to const Vector, Pressure") {
  quantities<Vector, Pressure> q{};
  quantities_ref<const Vector, Pressure> ref(q);
  REQUIRE(&ref[density] == &q[density]);
  REQUIRE(sizeof(ref) == 2 * sizeof(void*));
  REQUIRE(q[density] == 0);
  REQUIRE(ref[density] == 0);
  q[density] = 24;
  REQUIRE(q[density] == 24);
  REQUIRE(ref[density] == 24);
}

TEST_CASE("Can be assigned to") {
  quantities<Vector, Pressure> a{};
  quantities<Vector, Pressure> b({{}, 42});
  quantities_ref<Vector, Pressure> ref(a);
  REQUIRE(a != b);
  ref = b;
  REQUIRE(a == b);
}
