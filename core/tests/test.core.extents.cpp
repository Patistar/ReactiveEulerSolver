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

#include "fub/extents.hpp"

int main() {
  using fub::dynamic_extent;
  using fub::extents;

  constexpr extents<2, 2> e_static;
  static_assert(size(e_static) == 4, "Static Size is broken.");
  static_assert(e_static.extent(0) == 2, "Static Access is broken.");
  static_assert(e_static.extent(1) == 2, "Static Access is broken.");

  constexpr extents<2, dynamic_extent> e_mixed(2);
  static_assert(e_static == e_mixed, "");
  static_assert(e_mixed.extent(0) == 2, "");
  static_assert(e_mixed.extent(1) == 2, "");

  constexpr extents<dynamic_extent, dynamic_extent> e_dynamic(2, 2);
  static_assert(e_static == e_dynamic, "");
  static_assert(e_mixed == e_dynamic, "");
  static_assert(e_dynamic.extent(0) == 2, "");
  static_assert(e_dynamic.extent(1) == 2, "");

  constexpr std::array<std::ptrdiff_t, 2> array = as_array(e_mixed);
  static_assert(array[0] == 2, "");
  static_assert(array[1] == 2, "");
  constexpr std::array<std::ptrdiff_t, 1> dynamic =
      get_dynamic_extents(e_mixed);
  static_assert(dynamic[0] == 2, "");

  constexpr extents<2, dynamic_extent> grow_1_e_mixed = fub::grow<1>(e_mixed);
  static_assert(grow_1_e_mixed == extents<2, 3>(), "");

  constexpr extents<3, dynamic_extent> grow_0_e_mixed = fub::grow<0>(e_mixed);
  static_assert(grow_0_e_mixed == extents<3, 2>(), "");

  constexpr extents<dynamic_extent, dynamic_extent> grow_1_e_dynamic =
      fub::grow<1>(e_dynamic);
  static_assert(grow_1_e_dynamic == extents<2, 3>(), "");

  constexpr extents<dynamic_extent, dynamic_extent> grow_0_e_dynamic =
      fub::grow<0>(e_dynamic);
  static_assert(grow_0_e_dynamic == extents<3, 2>(), "");

  constexpr extents<2, 3> grow_1_e_static = fub::grow<1>(e_static);
  static_assert(grow_1_e_static == extents<2, 3>(), "");

  constexpr extents<3, 2> grow_0_e_static = fub::grow<0>(e_static);
  static_assert(grow_0_e_static == extents<3, 2>(), "");

  constexpr extents<4, dynamic_extent> replace_0 = fub::replace_extent<0, 4>(e_mixed);
  static_assert(replace_0 == extents<4, 2>(), "");

  constexpr extents<2, 4> replace_1 = fub::replace_extent<1, 4>(e_mixed);
  static_assert(replace_1 == extents<2, 4>(), "");
}
