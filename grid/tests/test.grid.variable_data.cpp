#include "fub/variable_data.hpp"

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

struct Density : fub::scalar_variable {};
struct Temperature : fub::scalar_variable {};
struct Variables : fub::variable_list<Density, Temperature> {};

constexpr auto density = fub::tag<Density>;
constexpr auto temperature = fub::tag<Temperature>;

TEST_CASE("Construct an interval of statically sized variable data") {
  fub::basic_variable_data<Variables, fub::mdspan<double, 10>> data{};
  SECTION("basic properties") {
    REQUIRE(data.size() == 20);
    REQUIRE(data.get_extents().extent(0) == 10);
  }
  SECTION("indexable through variables and index mapping") {
    fub::for_each_index(data.get_mapping(), [&](std::ptrdiff_t i) {
      data[density](i) = i;
      data[temperature](i) = 10 + i;
      REQUIRE(data(i)[density] == data[density](i));
      REQUIRE(data(i)[temperature] == data[temperature](i));
    });
    fub::for_each_index(data.get_mapping(), [&](std::ptrdiff_t i) {
      REQUIRE(data[density](i) == i);
      REQUIRE(data[temperature](i) == 10 + i);
    });
    auto variables = data.get_variable_list();
    fub::for_each(variables, [&](auto variable) {
      for (int i = 0; i < data.get_extents().extent(0); ++i) {
        REQUIRE(data[variable](i) == 10 * (*variables.index(variable)) + i);
      }
    });
  }
}