// Copyright (c) 2018 Maikel Nadolski
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "fub/span.hpp"

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

TEST_CASE("sizeof spans are small") {
  REQUIRE(sizeof(fub::span<int>) == sizeof(int*) + sizeof(std::ptrdiff_t));
  REQUIRE(sizeof(fub::span<int, 1>) == sizeof(int*));
}

TEST_CASE("Create an empty span.")
{
  fub::span<int> s;
  REQUIRE(s.data() == nullptr);
  REQUIRE(s.size() == 0);
}

TEST_CASE("Is contextually convertible to bool.") {
  SECTION("Dynamically-sized span") {
    // Default construction creates an empty span.
    fub::span<int> s;
    REQUIRE(!s);

    // After assignments the conversion returns true.
    int array[1];
    s = fub::span<int>{array};
    REQUIRE(s);
  }

  SECTION("Fixed-sized span") {
    // Fixed-sized spans can not be empty.
    int array[1];
    fub::span<int, 1> s{array};
    REQUIRE(s);
  }
}

