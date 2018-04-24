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

const int* take_const_span(fub::span<const int> x) {
  return x.data();
}

int* take_span(fub::span<int> x) {
  return x.data();
}

const int* take_const_span_1(fub::span<const int, 1> x) {
  return x.data();
}

int* take_span_1(fub::span<int, 1> x) {
  return x.data();
}

TEST_CASE("Create some empty spans")
{
  fub::span<int> s;
  fub::span<int, 1> t;

  REQUIRE(!s);
  REQUIRE(!t);
  REQUIRE(!take_const_span(s));
  REQUIRE(!take_span(s));
  REQUIRE(!take_const_span(t));
  REQUIRE(!take_span(t));
  REQUIRE(!take_const_span_1(t));
  REQUIRE(!take_span_1(t));
}

