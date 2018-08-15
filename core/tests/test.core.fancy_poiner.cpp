#include "fub/address_native.hpp"
#include "fub/fancy_pointer.hpp"

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

template <typename T> void test_nullable_concept_requirements() {
  using pointer = fub::fancy_pointer<T, fub::address_native>;
  REQUIRE(std::is_default_constructible<pointer>::value);
  REQUIRE(std::is_copy_constructible<pointer>::value);
  REQUIRE(std::is_copy_assignable<pointer>::value);
  REQUIRE(std::is_destructible<pointer>::value);
}

TEST_CASE("fancy_pointer<T, address_native> fulfills NullablePointer") {
  test_nullable_concept_requirements<int>();
  test_nullable_concept_requirements<void>();
  struct A;
  test_nullable_concept_requirements<A>();

  fub::fancy_pointer<int, fub::address_native> ptr = nullptr;
  REQUIRE(fub::fancy_pointer<int, fub::address_native>() == nullptr);
  REQUIRE(fub::fancy_pointer<int, fub::address_native>() == ptr);
  REQUIRE(ptr == fub::fancy_pointer<int, fub::address_native>());
}

TEST_CASE("fancy_pointer<T, address_native> gets correctly rebound.") {
  struct A;
  using rebound = typename std::pointer_traits<
      fub::fancy_pointer<int, fub::address_native>>::template rebind<A>;
  REQUIRE(
      std::is_same<rebound, fub::fancy_pointer<A, fub::address_native>>::value);
}

TEST_CASE("Point to some object") {
  int x = 42;
  fub::fancy_pointer<int, fub::address_native> p(&x);
  REQUIRE(*p == 42);
  *p = 24;
  REQUIRE(x == 24);
  int* xptr = std::addressof(*p);
  REQUIRE(xptr == &x);
}

TEST_CASE("Get distance of objects in an array") {
  int x[100];
  fub::fancy_pointer<int, fub::address_native> origin(&x[0]);
}