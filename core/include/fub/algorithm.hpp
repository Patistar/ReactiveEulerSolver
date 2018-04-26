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

#ifndef FUB_CORE_ALGORITHM_HPP
#define FUB_CORE_ALGORITHM_HPP

#include "fub/functional.hpp"
#include "fub/type_traits.hpp"
#include "fub/utility.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include <range/v3/range_concepts.hpp>
#include <range/v3/utility/invoke.hpp>
#include <type_traits>

namespace fub {
#ifdef FUB_CORE_USE_STD_CLAMP
using std::clamp;
#else
template <class T, class Compare = std::less<T>>
constexpr std::enable_if_t<std::is_arithmetic<T>::value, const T&>
clamp(const T& v, const T& lo, const T& hi, Compare comp = Compare()) {
  return assert(!comp(hi, lo)), comp(v, lo) ? lo : comp(hi, v) ? hi : v;
}
#endif

/// @brief Constexpr version of std::accumulate.
template <typename I, typename S, typename T, typename BinaryOp = std::plus<>>
constexpr std::enable_if_t<
    ranges::InputIterator<I>() && ranges::Sentinel<S, I>(), T>
accumulate(I first, S last, T init, BinaryOp binary_op = BinaryOp()) {
  while (first != last) {
    init = ranges::invoke(binary_op, init, *first++);
  }
  return init;
}

/// @brief Ranged constexpr version of accumulate.
template <typename R, typename T, typename BinaryOp = std::plus<>>
constexpr std::enable_if_t<ranges::Range<R>::value, T>
accumulate(R&& rng, T&& init, BinaryOp&& binary_op = BinaryOp()) {
  return fub::accumulate(ranges::begin(rng), ranges::end(rng),
                         std::forward<T>(init),
                         std::forward<BinaryOp>(binary_op));
}

template <typename L>
constexpr std::ptrdiff_t count(std::initializer_list<L>&& list,
                               const nodeduce_t<L>& needle) {
  auto first = list.begin();
  auto last = list.end();
  std::ptrdiff_t counter{0};
  while (first != last) {
    if (*first == needle) {
      counter += 1;
    }
    ++first;
  }
  return counter;
}

template <typename I, typename S, typename T, typename Pred = std::equal_to<>>
constexpr std::enable_if_t<
    conjunction<ranges::InputIterator<I>, ranges::Sentinel<S, I>>::value,
    std::ptrdiff_t>
count(I first, S last, const T& needle, Pred pred = Pred()) {
  std::ptrdiff_t counter{0};
  while (first != last) {
    if (ranges::invoke(pred, *first, needle)) {
      ++counter;
    }
    ++first;
  }
  return counter;
}

template <typename R, typename T, typename Pred = std::equal_to<>>
constexpr std::enable_if_t<ranges::InputRange<R>::value, std::ptrdiff_t>
count(R&& rng, const T& needle, Pred pred = Pred()) {
  return count(ranges::begin(rng), ranges::end(rng), needle, std::move(pred));
}

template <typename I, typename S, typename J, typename T,
          typename BinaryOp1 = std::plus<>,
          typename BinaryOp2 = std::multiplies<>>
constexpr std::enable_if_t<ranges::InputIterator<I>() &&
                               ranges::Sentinel<S, I>() &&
                               ranges::InputIterator<J>(),
                           T>
transform_reduce(I left, S last, J right, T initial,
                 BinaryOp1 binary_op1 = BinaryOp1(),
                 BinaryOp2 binary_op2 = BinaryOp2()) {
  while (left != last) {
    initial = ranges::invoke(binary_op1, initial,
                             fub::invoke(binary_op2, *left, *right));
    ++left;
    ++right;
  }
  return initial;
}

template <typename R, typename S, typename T, typename BinaryOp1 = std::plus<>,
          typename BinaryOp2 = std::multiplies<>>
constexpr std::enable_if_t<ranges::InputRange<R>() && ranges::InputRange<S>(),
                           T>
transform_reduce(R&& left, S&& right, T initial,
                 BinaryOp1 binary_op1 = BinaryOp1(),
                 BinaryOp2 binary_op2 = BinaryOp2()) {
  return fub::transform_reduce(ranges::begin(left), ranges::end(left),
                               ranges::begin(right), std::move(initial),
                               std::move(binary_op1), std::move(binary_op2));
}

template <typename A, typename I>
constexpr std::enable_if_t<std::is_arithmetic<A>::value, bool>
almost_equal(const A& x, const nodeduce_t<A>& y, const I& ulp) noexcept {
  constexpr A eps = std::numeric_limits<A>::epsilon();
  constexpr A min = std::numeric_limits<A>::min();
  const A diff = std::abs(x - y);
  const A sum = std::abs(x + y);
  return diff <= eps * sum * ulp || diff < min;
}

} // namespace fub

#endif // !ALGORITHM_HPP
