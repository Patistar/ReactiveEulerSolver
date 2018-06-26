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

#include "fub/octree.hpp"

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include <array>

TEST_CASE("Parent and select child are working together") {
  fub::octant<3> root{};
  fub::octant<3> child_1 = root.child(0);
  fub::octant<3> child_2 = root.child(6);
  fub::octant<3> child_3 = child_2.child(5);

  SECTION("Test with parent") {
    REQUIRE(child_1 != child_2);
    REQUIRE(child_1.parent() == child_2.parent());
    REQUIRE(root == child_2.parent());
  }
  SECTION("test inter child relation") {
    REQUIRE(child_2 != child_3);
    REQUIRE(child_2 == child_3.parent());
  }
}

template <int Rank> void test_max_depth() {
  constexpr int max_depth = fub::octant<Rank>::max_depth;
  std::array<fub::octant<Rank>, max_depth + 1> children{};
  for (int depth = 1; static_cast<std::size_t>(depth) < children.size();
       depth++) {
    children[depth] = children[depth - 1].child(1);
    REQUIRE(children[depth] != children[depth - 1]);
    REQUIRE(children[depth].depth() == depth);
    REQUIRE(children[depth].parent() == children[depth - 1]);
  }
}

TEST_CASE("max_depth is maximal depth") {
  SECTION("three dimensions") { test_max_depth<3>(); }
  SECTION("two dimensions") { test_max_depth<2>(); }
  SECTION("one dimensions") { test_max_depth<1>(); }
}

TEST_CASE("Construction with coordinates works with coordinate function") {
  SECTION("ones") {
    fub::octant<3> octant{4, {1, 1, 1}};
    REQUIRE(octant.depth() == 4);
    REQUIRE(fub::coordinate<0>(octant) == 1);
    REQUIRE(fub::coordinate<1>(octant) == 1);
    REQUIRE(fub::coordinate<2>(octant) == 1);
  }
  SECTION("ones coordinates") {
    std::array<std::uint64_t, 3> coords{1, 1, 1};
    fub::octant<3> octant{4, coords};
    REQUIRE(fub::coordinates(octant) == coords);

    coords = std::array<std::uint64_t, 3>{14, 7, 5};
    octant = fub::octant<3>{4, coords};
    REQUIRE(fub::coordinates(octant) == coords);
  }
  SECTION("some") {
    fub::octant<3> octant{4, {14, 7, 5}};
    REQUIRE(octant.depth() == 4);
    REQUIRE(fub::coordinate<0>(octant) == 14);
    REQUIRE(fub::coordinate<1>(octant) == 7);
    REQUIRE(fub::coordinate<2>(octant) == 5);
  }
}

TEST_CASE("face Neighbor") {
  std::array<std::uint64_t, 3> coordinates{13, 7, 5};
  fub::octant<3> octant{4, coordinates};
  SECTION("x axis") {
    const int x0 = 14;
    fub::face face{fub::axis::x, fub::direction::left};
    while (coordinates[0] != 0) {
      fub::optional<fub::octant<3>> nb = fub::face_neighbor(octant, face);
      coordinates[0] -= 1;
      REQUIRE(bool(nb));
      REQUIRE(fub::coordinates(*nb) == coordinates);
      octant = *nb;
    }
  }
}

constexpr std::uint64_t max_coordinate(int depth) noexcept {
  return depth ? (std::uint64_t(1) << depth) - 1 : 0;
}

void test_descendant(const fub::octant<3>& octant, int depth) {
  const std::uint64_t max = max_coordinate(depth - octant.depth());
  const auto coords = coordinates(octant);
  fub::octant<3> expected_lower{depth, coords};
  std::remove_const_t<decltype(coords)> upper_coords;
  std::transform(coords.begin(), coords.end(), upper_coords.begin(),
                 [max](auto x) { return x + max; });
  fub::octant<3> expected_upper{depth, upper_coords};
  fub::octant<3> lower = octant.lower_descendant(depth);
  fub::octant<3> upper = octant.upper_descendant(depth);
  REQUIRE(lower == expected_lower);
  REQUIRE(coordinates(upper) == coordinates(expected_upper));
}

TEST_CASE("descendants") {
  SECTION("construct bounds per hand") {
    for (int depth = 0; depth < 10; ++depth) {
      test_descendant(fub::octant<3>{}, depth);
    }
  }

  SECTION("count elements") {
    const auto first = fub::octant<3>{}.lower_descendant(3);
    const auto last = fub::octant<3>{}.upper_descendant(3);
    int counter{1};
    auto oct{first};
    while (oct != last) {
      ++counter;
      oct = oct.next();
    }
    REQUIRE(counter == (1 << 3) * (1 << 3) * (1 << 3));
  }
}

