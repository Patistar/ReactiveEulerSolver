#include "fub/variable_accessor.hpp"

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

struct Density {};
struct Pressure {};
struct Conservative : fub::basic_variable_list<Density, Pressure> {};

static_assert(fub::variable_traits<Conservative>::size() == 2, "Size is not two!");

TEST_CASE("to_pointer_storage usage") {
  using Accessor = fub::variable_accessor<Conservative>;
  Accessor accessor;
  Accessor::value_storage storage;
  Accessor::pointer_storage pointer = accessor.to_pointer_storage(storage);
  REQUIRE(std::get<0>(pointer) == std::get<0>(storage).data());
}