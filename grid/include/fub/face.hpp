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

#ifndef FUB_FACE_HPP
#define FUB_FACE_HPP

#include <type_traits>

namespace fub {

enum class direction { left = -1, right = +1 };

enum class axis { x, y, z };

template <axis A> using axis_constant = std::integral_constant<axis, A>;

template <axis A> static constexpr axis_constant<A> axis_c{};

struct face {
  axis dimension;
  direction side;
};

constexpr int as_int(axis dim) noexcept { return static_cast<int>(dim); }
constexpr int sign(direction side) noexcept { return static_cast<int>(side); }

template <typename Archive>
void save(Archive& archive, const face& face, unsigned /* version */) {
  int dim = static_cast<int>(face.dimension);
  int side = static_cast<int>(face.side);
  archive << dim;
  archive << side;
}

template <typename Archive>
void load(Archive& archive, face& face, unsigned /* version */) {
  int dim;
  int side;
  archive << dim;
  archive << side;
  face.dimension = static_cast<axis>(dim);
  face.side = static_cast<direction>(side);
}

} // namespace fub

#endif // !FACE_HPP
