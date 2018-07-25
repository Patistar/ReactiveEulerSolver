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

#include "fub/simd_span.hpp"
#include <array>

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

TEST_CASE("Create a simd span mutable") {
  using simd_type = fub::simd<double>;
  alignas(fub::memory_alignment_v<simd_type>) std::array<double, 8> array{{1, 2, 3, 4, 5, 6, 7, 8}};
  fub::simd_span<double> span(array.data(), 4);
  REQUIRE(span.size() == array.size() / fub::simd_size_v<double>);
  REQUIRE(array[0] == 1);
  span[0] = simd_type(42);
  for (int i = 0; i < fub::simd_size_v<double>; ++i) {
    REQUIRE(array[i] == 42);
  }
  REQUIRE(std::is_assignable<simd_type, decltype(span[0])>::value);
}

TEST_CASE("Create a simd span to const") {
  using simd_type = fub::simd<double>;
  alignas(fub::memory_alignment_v<simd_type>) std::array<double, 8> array{{1, 2, 3, 4, 5, 6, 7, 8}};
  fub::simd_span<const double> span(array.data(), 4);
  REQUIRE(span.size() == array.size() / fub::simd_size_v<double>);
  REQUIRE(array[0] == 1);
  REQUIRE(!std::is_assignable<decltype(span[0]), simd_type>::value);
  const fub::simd<double> s = span[0];
  for (int i = 0; i < s.size(); ++i) {
    REQUIRE(s[i] == array[i]);
  }
}