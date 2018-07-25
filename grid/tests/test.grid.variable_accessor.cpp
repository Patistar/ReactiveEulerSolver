#include "fub/accessor_simd.hpp"
#include "fub/variable_accessor.hpp"

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include <boost/hana.hpp>

struct Density {};
struct Pressure {};
template <int Dim> struct Velocity {};
struct Temperature {};
struct Primitive : fub::variable_union<Temperature, Pressure> {};

using boost::hana::nothing;
using boost::hana::type_c;
using fub::types;

static_assert(fub::variable_traits<Primitive>::size() == 2, "Size is not two!");

template <typename T>
using finder_t = fub::detail::variable_accessor_find_index<T>;

static_assert(finder_t<Density>::find_union(type_c<Density>) != nothing,
              "find_union is broken.");

static_assert(finder_t<const Density>::find_union(type_c<Density>) != nothing,
              "find_union is broken.");

static_assert(finder_t<Primitive>::find_union(type_c<Temperature>) != nothing,
              "find_union is broken.");

static_assert(finder_t<const Primitive>::find_union(type_c<Temperature>) !=
                  nothing,
              "find_union is broken.");

static_assert(finder_t<Primitive>::is_union(type_c<Primitive>),
              "is_union is broken.");

static_assert(finder_t<const Primitive>::is_union(type_c<Primitive>),
              "is_union is broken.");

static_assert(finder_t<Primitive>::is_union(type_c<const Primitive>),
              "is_union is broken.");

TEST_CASE("test for Density typedefs") {
  using ConstAccessor = fub::variable_accessor<Density>;
  REQUIRE(type_c<ConstAccessor::element_type<Density>> == type_c<double>);
  REQUIRE(type_c<ConstAccessor::element_type<const Density>> ==
          type_c<const double>);
  REQUIRE(type_c<ConstAccessor::value_type<const Density>> == type_c<double>);
  REQUIRE(type_c<ConstAccessor::value_type<Density>> == type_c<double>);
}

TEST_CASE("test for const Density typedefs") {
  using ConstAccessor = fub::variable_accessor<const Density>;
  REQUIRE(type_c<ConstAccessor::element_type<Density>> == type_c<const double>);
  REQUIRE(type_c<ConstAccessor::element_type<const Density>> ==
          type_c<const double>);
  REQUIRE(type_c<ConstAccessor::value_type<const Density>> == type_c<double>);
  REQUIRE(type_c<ConstAccessor::value_type<Density>> == type_c<double>);
}

TEST_CASE("test for types<Density> typedefs") {
  using ConstAccessor = fub::variable_accessor<types<Density>>;
  REQUIRE(type_c<ConstAccessor::element_type<Density>> == type_c<double>);
  REQUIRE(type_c<ConstAccessor::element_type<const Density>> ==
          type_c<const double>);
  REQUIRE(type_c<ConstAccessor::value_type<const Density>> == type_c<double>);
  REQUIRE(type_c<ConstAccessor::value_type<Density>> == type_c<double>);
}

TEST_CASE("test for types<const Density> typedefs") {
  using ConstAccessor = fub::variable_accessor<types<const Density>>;
  REQUIRE(type_c<ConstAccessor::element_type<Density>> == type_c<const double>);
  REQUIRE(type_c<ConstAccessor::element_type<const Density>> ==
          type_c<const double>);
  REQUIRE(type_c<ConstAccessor::value_type<const Density>> == type_c<double>);
  REQUIRE(type_c<ConstAccessor::value_type<Density>> == type_c<double>);
}

TEST_CASE("test for Primitive typedefs") {
  using ConstAccessor = fub::variable_accessor<Primitive>;
  REQUIRE(type_c<ConstAccessor::element_type<Temperature>> == type_c<double>);
  REQUIRE(type_c<ConstAccessor::element_type<const Temperature>> ==
          type_c<const double>);
  REQUIRE(type_c<ConstAccessor::value_type<const Temperature>> ==
          type_c<double>);
  REQUIRE(type_c<ConstAccessor::value_type<Temperature>> == type_c<double>);
}

TEST_CASE("test for const Primitive typedefs") {
  using ConstAccessor = fub::variable_accessor<const Primitive>;
  REQUIRE(type_c<ConstAccessor::element_type<Temperature>> ==
          type_c<const double>);
  REQUIRE(type_c<ConstAccessor::element_type<const Temperature>> ==
          type_c<const double>);
  REQUIRE(type_c<ConstAccessor::value_type<const Temperature>> ==
          type_c<double>);
  REQUIRE(type_c<ConstAccessor::value_type<Temperature>> == type_c<double>);
}

TEST_CASE("test for types<const Primitive, Density> typedefs") {
  using ConstAccessor = fub::variable_accessor<types<const Primitive, Density>>;
  REQUIRE(type_c<ConstAccessor::element_type<Temperature>> ==
          type_c<const double>);
  REQUIRE(type_c<ConstAccessor::element_type<const Temperature>> ==
          type_c<const double>);
  REQUIRE(type_c<ConstAccessor::value_type<const Temperature>> ==
          type_c<double>);
  REQUIRE(type_c<ConstAccessor::value_type<Temperature>> == type_c<double>);
  REQUIRE(type_c<ConstAccessor::element_type<Density>> == type_c<double>);
  REQUIRE(type_c<ConstAccessor::element_type<const Density>> ==
          type_c<const double>);
  REQUIRE(type_c<ConstAccessor::value_type<const Density>> == type_c<double>);
  REQUIRE(type_c<ConstAccessor::value_type<Density>> == type_c<double>);
  REQUIRE(type_c<ConstAccessor::span<const Density>> == type_c<fub::span<const double>>);
  REQUIRE(type_c<ConstAccessor::span<Density>> == type_c<fub::span<double>>);
  REQUIRE(type_c<ConstAccessor::span<const Temperature>> == type_c<fub::span<const double>>);
  REQUIRE(type_c<ConstAccessor::span<Temperature>> == type_c<fub::span<const double>>);
}

TEST_CASE("storage for single variable") {
  using Accessor = fub::variable_accessor<Density>;
  Accessor accessor;
  Accessor::value_storage storage{42.};
  Accessor::pointer_storage pointer = accessor.to_pointer_storage(storage);
  REQUIRE(pointer == &storage);
  auto rho = accessor.access(type_c<Density>, storage);
  REQUIRE(rho == storage);
  auto vars = accessor.get_accessible_variables();
  static_assert(boost::hana::size(vars) == boost::hana::size_c<1>, "");
}

TEST_CASE("storage for union") {
  using Accessor = fub::variable_accessor<Primitive>;
  Accessor accessor;
  Accessor::value_storage storage{42., 24.};
  Accessor::pointer_storage pointer = accessor.to_pointer_storage(storage);
  REQUIRE(pointer == &storage[0]);
  auto T = accessor.access(type_c<Temperature>, fub::as_const(storage), 0);
  REQUIRE(T == storage[0]);
  auto p = accessor.access(type_c<Pressure>, storage);
  REQUIRE(p == storage[1]);
  auto vars = accessor.get_accessible_variables();
  static_assert(boost::hana::size(vars) == boost::hana::size_c<2>, "");
}

TEST_CASE("storage for multiple unions") {
  using Accessor = fub::variable_accessor<fub::types<Primitive, Density>>;
  constexpr Accessor accessor;
  constexpr Accessor::value_storage storage{{42., 24.}, {12.}};
  auto vars = accessor.get_accessible_variables();
  // static_assert(boost::hana::size(vars) == boost::hana::size_c<3>, "");
  const double& T =
      accessor.access(type_c<const Temperature>, storage[fub::int_c<0>], 0);
  auto rho = accessor.access(type_c<Density>, storage);
  auto p = accessor.access(type_c<const Pressure>, storage[fub::int_c<0>], 0);
  REQUIRE(T == 42.0);
  REQUIRE(p == 24.0);
  REQUIRE(rho == 12.0);
  // Accessor::pointer_storage pointer = accessor.to_pointer_storage(storage);
}

TEST_CASE("storage for simd variables") {
  using Accessor = fub::variable_accessor<
      fub::types<Primitive>,
      fub::accessor_simd_aligned<double, fub::simd_abi::sse>>;
  Accessor accessor;
  Accessor::value_storage storage{{42., 24.}};
  {
    auto temp = accessor.access(type_c<Temperature>, storage);
    auto index = accessor.index(type_c<Temperature>);
    REQUIRE(all_of(temp == storage[index][0]));
  }
  {
    auto p = accessor.access(type_c<Pressure>, storage);
    auto index = accessor.index(type_c<Pressure>);
    REQUIRE(all_of(p == storage[index][1]));
  }
  auto vars = accessor.get_accessible_variables();
  static_assert(boost::hana::size(vars) == boost::hana::size_c<2>, "");
}