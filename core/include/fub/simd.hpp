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

#ifndef FUB_CORE_SIMD_HPP
#define FUB_CORE_SIMD_HPP

#include <Vc/vector.h>

#include "fub/tuple.hpp"
#include "fub/utility.hpp"

#include <array>
#include <cassert>
#include <limits>
#include <ostream>

namespace fub {
inline namespace v1 {

template <typename T>
Vc::Vector<T> clamp(const Vc::Vector<T>& v, const nodeduce_t<Vc::Vector<T>>& lo,
                   const nodeduce_t<Vc::Vector<T>>& hi) {
  assert(all_of(lo < hi));
  Vc::Vector<T> x{v};
  where(x < lo, x) = lo;
  where(hi < x, x) = hi;
  return x;
}

inline bool all_of(bool mask) noexcept { return mask; }
inline bool any_of(bool mask) noexcept { return mask; }

template <typename T>
struct where_expression {
  where_expression() = delete;
  where_expression(const where_expression&) = delete;
  where_expression(where_expression&&) = default;

  void operator=(const T& other) const {
    if (m_mask) {
      *m_data = other;
    }
  }

  bool m_mask;
  T* m_data;
};

template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>{}>>
where_expression<T> where(bool mask, T& data) {
  return where_expression<T>{mask, &data};
}

} // namespace v1
} // namespace fub

#endif
