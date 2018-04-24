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

#include "fub/patch.hpp"

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

struct Density {
  using value_type = double;
};

struct Velocity {
  using value_type = int;
};

struct Pressure {
  using value_type = double;
};

using fub::automatic_storage_descriptor;
using fub::dynamic_storage_descriptor;
using fub::extents;
using fub::patch;
using fub::pmr_storage_descriptor;
using fub::mdspan;

TEST_CASE("2D patch with simple variables") {
  patch<std::tuple<Density, Velocity, Pressure>, extents<16, 16>> block{};

  auto density = block.get<Density>();
  auto velocity = block.get<Velocity>();
  auto pressure = block.get<Pressure>();
  REQUIRE(
      std::is_same<mdspan<double, extents<16, 16>>, decltype(density)>::value);
  REQUIRE(density.size() == 16 * 16);
  REQUIRE(velocity.size() == 16 * 16);
  REQUIRE(pressure.size() == 16 * 16);
}

TEST_CASE("Using the automatic allocator") {
  using Vars = std::tuple<Density, Velocity, Pressure>;
  using Storage = automatic_storage_descriptor;
  patch<Vars, extents<16, 16>, Storage> block{};
  auto density = block.get<Density>();
  REQUIRE(
      std::is_same<mdspan<double, extents<16, 16>>, decltype(density)>::value);
  REQUIRE(density.size() == 16 * 16);
}

TEST_CASE("Using the dynamic allocator") {
  using Vars = std::tuple<Density, Velocity, Pressure>;
  using Storage = dynamic_storage_descriptor<>;
  patch<Vars, extents<16, 16>, Storage> block{};
  auto density = block.get<Density>();
  REQUIRE(
      std::is_same<mdspan<double, extents<16, 16>>, decltype(density)>::value);
  REQUIRE(density.size() == 16 * 16);
}

struct MyMemoryResource : boost::container::pmr::memory_resource {
  int num_allocates{0};
  int num_deallocates{0};

  void* do_allocate(std::size_t size, std::size_t alignment) override {
    boost::container::pmr::memory_resource* upstream =
        boost::container::pmr::get_default_resource();
    ++num_allocates;
    return upstream->allocate(size, alignment);
  }

  void do_deallocate(void* ptr, std::size_t size,
                     std::size_t alignment) override {
    boost::container::pmr::memory_resource* upstream =
        boost::container::pmr::get_default_resource();
    ++num_deallocates;
    return upstream->deallocate(ptr, size, alignment);
  }

  bool do_is_equal(const memory_resource& other) const noexcept override {
    return dynamic_cast<const MyMemoryResource*>(&other);
  }
};

TEST_CASE("Using the pmr allocator") {
  using Vars = std::tuple<Density, Velocity, Pressure>;
  using Storage = pmr_storage_descriptor;
  MyMemoryResource resource{};
  Storage descriptor{{&resource}};
  {
    patch<Vars, extents<16, 16>, Storage> block{descriptor};
    REQUIRE(resource.num_allocates == 3);
    REQUIRE(resource.num_deallocates == 0);
    auto density = block.get<Density>();
    REQUIRE(
        std::is_same<mdspan<double, extents<16, 16>>, decltype(density)>::value);
    REQUIRE(density.size() == 16 * 16);
  }
  REQUIRE(resource.num_allocates == 3);
  REQUIRE(resource.num_deallocates == 3);
}
