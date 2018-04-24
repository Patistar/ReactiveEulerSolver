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

#ifndef FUB_TUPLE_HPP
#define FUB_TUPLE_HPP

#include "fub/algorithm.hpp"
#include "fub/functional.hpp"

#include <range/v3/utility/tuple_algorithm.hpp>
#include <tuple>

namespace fub {
template <typename... Tuples>
using tuple_cat_t = decltype(std::tuple_cat(std::declval<Tuples>()...));
}

#ifdef FUB_CORE_USE_STD_APPLY
namespace fub {
using std::apply;
}
#else
namespace fub {
class tuple_apply_fn {
  template <typename F, typename T, std::size_t... Is>
  constexpr decltype(auto) impl(F&& fun, T&& tuple,
                                std::index_sequence<Is...>) {
    return fub::invoke(std::forward<F>(fun), std::get<Is>(tuple)...);
  }

public:
  template <typename F, typename T>
  constexpr decltype(auto) operator()(F&& fun, T&& tuple) {
    return impl(
        std::forward<F>(fun), std::forward<T>(tuple),
        std::make_index_sequence<std::tuple_size<std::decay_t<T>>::value>());
  }
};

template <typename F, typename T>
constexpr decltype(auto) apply(F&& fun, T&& tuple) {
  return tuple_apply_fn{}(fun, tuple);
}
} // namespace fub
#endif
namespace fub {
template <typename Tuple, typename T, typename F>
T foldl(Tuple&& tuple, T&& init, F&& fun) noexcept {
  return ranges::tuple_foldl(std::forward<Tuple>(tuple), std::forward<T>(init),
                             std::forward<F>(fun));
}

template <typename F, typename... Ts>
constexpr F for_each_tuple_element(F f, Ts&&... ts) {
  (void)std::initializer_list<int>{
      ((void)invoke(std::ref(f), std::forward<Ts>(ts)), 42)...};
  return f;
}

template <typename F, typename... Ts>
constexpr F for_each_tuple_element(F f, std::tuple<Ts...>&& tuple) {
  fub::apply(
      [&](auto&&... ts) {
        (void)std::initializer_list<int>{
            ((void)invoke(std::ref(f), std::forward<decltype(ts)>(ts)), 42)...};
      },
      std::move(tuple));
  return f;
}

template <typename F, typename... Ts>
constexpr F for_each_tuple_element(F f, const std::tuple<Ts...>& tuple) {
  fub::apply(
      [&](auto&&... ts) {
        (void)std::initializer_list<int>{
            ((void)invoke(std::ref(f), std::forward<decltype(ts)>(ts)), 42)...};
      },
      std::move(tuple));
  return f;
}

template <typename F, typename... Ts>
constexpr F for_each_tuple_element(F f, std::tuple<Ts...>& tuple) {
  fub::apply(
      [&](auto&&... ts) {
        (void)std::initializer_list<int>{
            ((void)invoke(std::ref(f), std::forward<decltype(ts)>(ts)), 42)...};
      },
      tuple);
  return f;
}

#if __cpp_lib_make_from_tuple
using std::make_from_tuple;
#else
namespace algorithm_detail {
template <class T, class Tuple, std::size_t... I>
constexpr T make_from_tuple_impl(Tuple&& t, std::index_sequence<I...>) {
  return T(std::get<I>(std::forward<Tuple>(t))...);
}
} // namespace algorithm_detail

template <class T, class Tuple> constexpr T make_from_tuple(Tuple&& t) {
  return algorithm_detail::make_from_tuple_impl<T>(
      std::forward<Tuple>(t),
      std::make_index_sequence<
          std::tuple_size<std::remove_reference_t<Tuple>>::value>{});
}
#endif

} // namespace fub

#endif // !TUPLE_HPP
