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

#ifndef FUB_CORE_TUPLE_HPP
#define FUB_CORE_TUPLE_HPP

#include "fub/functional.hpp"

#include <tuple>

namespace fub {
inline namespace v1 {
#ifdef FUB_WITH_STD_APPLY
using std::apply;
using std::make_from_tuple;
#else
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
        std::make_index_sequence<std::tuple_size<remove_cvref_t<T>>::value>());
  }
};

template <typename F, typename T>
constexpr decltype(auto) apply(F&& fun, T&& tuple) {
  return tuple_apply_fn{}(fun, tuple);
}

template <typename T, typename Tuple>
constexpr T make_from_tuple(Tuple&& t) {
  return fub::apply([](auto&&... args) {
    return T(std::forward<decltype(args)>(args)...);
  }, t);
}
#endif
} // namespace v1
} // namespace fub

#endif // !FUB_CORE_TUPLE_HPP
