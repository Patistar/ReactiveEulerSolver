// Copyright (c) 2017-2018 Maikel Nadolski
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

#ifndef FUB_FACE_HPP
#define FUB_FACE_HPP

#include <type_traits>

namespace fub {
inline namespace v1 {
enum class direction { left = -1, right = +1 };

enum class axis { x, y, z };

template <axis A> using axis_constant = std::integral_constant<axis, A>;

template <axis A> static constexpr axis_constant<A> axis_c{};

struct face {
  face() = default;
  constexpr face(axis dim, direction dir) : dimension{dim}, side{dir} {}
  explicit face(int i) : dimension(axis(i / 2)), side(direction(i % 2)) {}

  axis dimension;
  direction side;

  operator int() const noexcept {
    const int dir = (this->side == direction::left) ? 0 : 1;
    const int dim = 2 * static_cast<int>(dimension);
    const int result = dim + dir;
    return result;
  }
};

constexpr int as_int(axis dim) noexcept { return static_cast<int>(dim); }
constexpr int as_int(face f) noexcept {
  const int side = (f.side == direction::left) ? 0 : 1;
  const int dim = 2 * as_int(f.dimension);
  return dim + side;
}
constexpr int sign(direction side) noexcept { return static_cast<int>(side); }

constexpr direction invert(direction side) noexcept {
  return side == direction::left ? direction::right : direction::left;
}

} // namespace v1
} // namespace fub

#endif // !FACE_HPP
