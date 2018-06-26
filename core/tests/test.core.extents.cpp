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
#include <cstdio>

int main() {
  using fub::dyn;
  {
    constexpr fub::extents<2, 2> e;
    static_assert(e.size() == 4);
    static_assert(e.get(0) == 2);
    static_assert(e.get(1) == 2);
  }

  {
    constexpr fub::extents<2, dyn> e(2);
    static_assert(e == fub::extents<2, 2>(), "Only values shall matter.");
  }

  {
    constexpr fub::extents<dyn, 2> e(2);
    static_assert(e == fub::extents<2, 2>(), "Only values shall matter.");
  }

  {
    constexpr fub::extents<dyn, dyn> e(2, 2);
    static_assert(e == fub::extents<2, 2>(), "Only values shall matter.");
  }

  {
    constexpr fub::extents<2, 2> e;
    constexpr auto e2 = fub::grow(e, fub::int_c<1>);
    static_assert(std::is_same<std::decay_t<decltype(e2)>, fub::extents<2, 3>>::value, "Types are wrong.");
    static_assert(e2 == fub::extents<2, 3>(), "Equality Operator failed, but types are correct.");
  }

  {
    constexpr fub::extents<dyn, 2> e(2);
    constexpr auto e2 = fub::grow(e, fub::int_c<1>);
    static_assert(std::is_same<std::decay_t<decltype(e2)>, fub::extents<dyn, 3>>::value, "Types are wrong.");
    static_assert(e2 == fub::extents<2, 3>(), "Equality Operator failed, but types are correct.");
  }

  {
    constexpr fub::extents<2, dyn> e(2);
    constexpr auto e2 = fub::grow(e, fub::int_c<1>);
    static_assert(
        std::is_same<std::decay_t<decltype(e2)>, fub::extents<2, dyn>>::value);
    static_assert(e2 == fub::extents<2, 3>());
  }
}
