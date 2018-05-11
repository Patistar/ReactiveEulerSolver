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

#include "fub/algorithm.hpp"
#include "fub/patch_view.hpp"

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include <iostream>

struct Density {
  using value_type = double;
};

struct Velocity {
  using value_type = int;
};
static constexpr Velocity velocity{};

struct Pressure {
  using value_type = double;
};

using Extents = fub::extents<16, fub::dyn>;
using Patch = fub::patch<std::tuple<Density, Velocity, Pressure>, Extents>;

template <typename... Vars> using View = fub::patch_view<Extents, Vars...>;

bool test_function(View<const Velocity> view) {
  auto u = view[velocity];
  REQUIRE(u.size() == 16 * 100);
  return true;
}

TEST_CASE("View selects its variables.") {
  Patch patch(Extents(100));
  REQUIRE(test_function(patch));
}

TEST_CASE("Rows") {
  Patch patch(Extents(100));
  auto view = fub::make_view(patch);
  REQUIRE(
      std::is_same<decltype(view), fub::patch_view<Extents, Density, Velocity,
                                                   Pressure>>::value);
  auto rows = view.rows();
  int counter = 0;
  SECTION("hand written for loop") {
    for (auto row : rows) {
      REQUIRE(
          std::is_same<decltype(row),
                       fub::row_view<16, Density, Velocity, Pressure>>::value);
      REQUIRE(row.size() == 16);
      ++counter;
      if (counter > 100) {
        break;
      }
    }
    REQUIRE(counter == 100);
  }
  SECTION("for_each_row") {
    fub::for_each_row(
        [&counter](auto row) {
          REQUIRE(
              std::is_same<decltype(row), fub::row_view<16, Density, Velocity,
                                                        Pressure>>::value);
          REQUIRE(row.size() == 16);
          ++counter;
          if (counter > 100) {
            throw;
          }
        },
        view);
    REQUIRE(counter == 100);
  }
}

TEST_CASE("subview of rows") {
  using Vars = std::tuple<Density, Velocity, Pressure>;
  static constexpr int size = 10;
  fub::patch<Vars, fub::extents<size, size>> patch;
  int counter = 1;
  fub::for_each_row(
      [&](const fub::view_row_t<Vars, size>& row) {
        for (auto x : row) {
          fub::for_each_tuple_element(
              [=](auto q) { x[q] = counter; },
              std::tuple<Density, Velocity, Pressure>{});
          ++counter;
        }
      },
      fub::make_view(patch));

  SECTION("take and drop") {
    counter = 1;
    fub::for_each_row(
        [&](const fub::view_row_t<Vars, size>& row) {
          for (int i = 0; i < size; ++i) {
            fub::for_each_tuple_element(
                [=](auto q) { REQUIRE(row[q](i) == counter + i); },
                std::tuple<Density, Velocity, Pressure>{});
          }
          auto left = fub::take<size - 1>(row);
          auto right = fub::drop<1>(row);
          REQUIRE(left.size() == size - 1);
          REQUIRE(right.size() == size - 1);
          for (int i = 0; i < size - 1; ++i) {
            fub::for_each_tuple_element(
                [=](auto q) {
                  REQUIRE(left[q](i) == counter + i);
                  REQUIRE(right[q](i) == counter + 1 + i);
                  REQUIRE(left(i)[q] == counter + i);
                  REQUIRE(right(i)[q] == counter + 1 + i);
                },
                std::tuple<Density, Velocity, Pressure>{});
          }
          counter += size;
        },
        fub::make_view(patch));
  }

  SECTION("join") {
    fub::for_each_row(
        [](const fub::view_row_t<Vars, size>& row) {
          auto left = fub::rdrop<1>(row);
          auto right = fub::drop<1>(row);
          auto joined = fub::join(left, right);
          REQUIRE(joined.size() == 2 * size - 2);
        },
        fub::make_view(patch));
  }
}

TEST_CASE("vector variables") {
  struct Thermo : fub::vector_variable<Density, Pressure> {};
  fub::patch<std::tuple<Thermo, Velocity>, fub::extents<100>> patch{};
  auto view = fub::make_view(patch);
  fub::for_each_index(view.extents(), [&](auto index) {
    view(index) = fub::quantities<Thermo, Velocity>{{1.0, 2.0}, static_cast<int>(index[0])};
  });

  static constexpr Density density{};
  static constexpr Pressure pressure{};
  static constexpr Velocity velocity{};
  fub::for_each_index(view.extents(), [&](auto index) {
    REQUIRE(view[density](index) == 1.0);
    REQUIRE(view[pressure](index) == 2.0);
    REQUIRE(view[velocity](index) == index[0]);
    REQUIRE(view(index)[density] == 1.0);
    REQUIRE(view(index)[pressure] == 2.0);
    REQUIRE(view(index)[velocity] == index[0]);
  });
}
