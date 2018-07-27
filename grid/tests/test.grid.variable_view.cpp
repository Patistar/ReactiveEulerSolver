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

#include "fub/variable_view.hpp"

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

struct Density : fub::scalar_variable {};
struct Pressure : fub::scalar_variable {};
struct Velocity : fub::vector_variable<2> {};
struct Species : fub::vector_variable<fub::dynamic_extent> {
  using fub::vector_variable<fub::dynamic_extent>::vector_variable;
};

struct PrimitiveVariables : fub::variable_list<Density, Velocity, Pressure> {};
struct Variables : fub::variable_list<Density, Velocity, Pressure, Species> {
  using fub::variable_list<Density, Velocity, Pressure, Species>::variable_list;
};

constexpr auto density = fub::tag<Density>;
constexpr auto pressure = fub::tag<Pressure>;
template <std::ptrdiff_t Dim> constexpr auto velocity = fub::tag<Velocity, Dim>;

template <typename Vs, std::ptrdiff_t... Es>
using patch_view = fub::variable_view<Vs, double, fub::extents<Es...>>;

TEST_CASE("Have a view on a contiguous array") {
  constexpr std::ptrdiff_t data_size =
      patch_view<PrimitiveVariables, 10, 10>::static_size();
  REQUIRE(data_size > 0);
  double data[data_size];
  patch_view<PrimitiveVariables, 10, 10> patch(PrimitiveVariables(), data);
  REQUIRE(patch.size() == data_size);
  fub::for_each_index(
      patch.get_mapping(), [&](std::ptrdiff_t i, std::ptrdiff_t j) {
        patch[density](i, j) = patch.get_mapping()(i, j);
        patch[velocity<0>](i, j) = 100 + patch.get_mapping()(i, j);
        patch[velocity<1>](i, j) = 200 + patch.get_mapping()(i, j);
        patch[pressure](i, j) = 300 + patch.get_mapping()(i, j);
      });
  fub::for_each_index(
      patch.get_mapping(), [&](std::ptrdiff_t i, std::ptrdiff_t j) {
        REQUIRE(patch[density](i, j) == patch.get_mapping()(i, j));
        REQUIRE(patch[velocity<0>](i, j) == 100 + patch.get_mapping()(i, j));
        REQUIRE(patch[velocity<1>](i, j) == 200 + patch.get_mapping()(i, j));
        REQUIRE(patch[pressure](i, j) == 300 + patch.get_mapping()(i, j));
      });
}

TEST_CASE("Have a view on a dynamically sized array") {
  using View = patch_view<Variables, 10, 10>;
  REQUIRE(View::static_size() == fub::dynamic_extent);
  constexpr Variables variables(Density(), Velocity(), Pressure(), Species(10));
  constexpr std::ptrdiff_t data_size = View::static_size(variables);
  std::vector<double> data(data_size);
  View patch(variables, data);
  REQUIRE(patch.size() == data_size);
  fub::for_each_index(
      patch.get_mapping(), [&](std::ptrdiff_t i, std::ptrdiff_t j) {
        std::ptrdiff_t offset = patch.get_mapping()(i, j);
        patch[density](i, j) = offset;
        patch[velocity<0>](i, j) = 100 + offset;
        patch[velocity<1>](i, j) = 200 + offset;
        patch[pressure](i, j) = 300 + offset;
        fub::for_each(Species(10), [&](auto species) {
          patch[species](i, j) = 400 + 100 * species.index() + offset;
        });
      });
  fub::for_each_index(
      patch.get_mapping(), [&](std::ptrdiff_t i, std::ptrdiff_t j) {
        std::ptrdiff_t offset = patch.get_mapping()(i, j);
        REQUIRE(patch[density](i, j) == offset);
        REQUIRE(patch[velocity<0>](i, j) == 100 + offset);
        REQUIRE(patch[velocity<1>](i, j) == 200 + offset);
        REQUIRE(patch[pressure](i, j) == 300 + offset);
        fub::for_each(Species(10), [&](auto species) {
          REQUIRE(patch[species](i, j) == 400 + 100 * species.index() + offset);
        });
      });
}