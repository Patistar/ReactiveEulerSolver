#include "fub/fancy_pointer.hpp"
#include "fub/native_address.hpp"

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

template <typename T> void test_nullable_concept_requirements() {
  using pointer = fub::fancy_pointer<T, fub::native_address>;
  REQUIRE(std::is_default_constructible<pointer>::value);
  REQUIRE(std::is_copy_constructible<pointer>::value);
  REQUIRE(std::is_copy_assignable<pointer>::value);
  REQUIRE(std::is_destructible<pointer>::value);
}

TEST_CASE("fancy_pointer<T, native_address> fulfills NullablePointer") {
  test_nullable_concept_requirements<int>();
  test_nullable_concept_requirements<void>();
  struct A;
  test_nullable_concept_requirements<A>();

  fub::fancy_pointer<int, fub::native_address> ptr = nullptr;
  REQUIRE(fub::fancy_pointer<int, fub::native_address>() == nullptr);
  REQUIRE(fub::fancy_pointer<int, fub::native_address>() == ptr);
  REQUIRE(ptr == fub::fancy_pointer<int, fub::native_address>());
}

TEST_CASE("fancy_pointer<T, native_address> gets correctly rebound.") {
  struct A;
  using rebound = typename std::pointer_traits<
      fub::fancy_pointer<int, fub::native_address>>::template rebind<A>;
  REQUIRE(
      std::is_same<rebound, fub::fancy_pointer<A, fub::native_address>>::value);
}


TEST_CASE("Point to some object") {
  int x = 42;
  fub::fancy_pointer<int, fub::native_address> p(&x);
  REQUIRE(*p == 42);
  *p = 24;
  REQUIRE(x == 24);
  int* xptr = std::addressof(*p);
  REQUIRE(xptr == &x);
}

TEST_CASE("Get distance of objects in an array") {
  int x[100];
  fub::fancy_pointer<int, fub::native_address> origin(&x[0]);
}