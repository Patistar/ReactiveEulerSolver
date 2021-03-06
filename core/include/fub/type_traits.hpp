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

/// @file This file adds basic type traits utilities which are not yet
/// implemented in all standard libraries.

#ifndef FUB_TYPE_TRAITS_HPP
#define FUB_TYPE_TRAITS_HPP

#include <meta/meta.hpp>
#include <range/v3/data.hpp>
#include <range/v3/size.hpp>
#include <range/v3/utility/functional.hpp>
#include <tuple>
#include <type_traits>

namespace fub {
////////////////////////////////////////////////////////////////////////////////
//                                                          [traits.is_detected]

template <class...> using void_t = void;

struct nonesuch {
  nonesuch() = delete;
  ~nonesuch() = delete;
  nonesuch(nonesuch const&) = delete;
  void operator=(nonesuch const&) = delete;
};

namespace traits_detail {
template <class Default, class AlwaysVoid, template <class...> class Op,
          class... Args>
struct detector {
  using value_t = std::false_type;
  using type = Default;
};

template <class Default, template <class...> class Op, class... Args>
struct detector<Default, void_t<Op<Args...>>, Op, Args...> {
  using value_t = std::true_type;
  using type = Op<Args...>;
};
} // namespace traits_detail

template <template <class...> class Op, class... Args>
using is_detected =
    typename traits_detail::detector<nonesuch, void, Op, Args...>::value_t;

template <template <class...> class Op, class... Args>
using detected_t =
    typename traits_detail::detector<nonesuch, void, Op, Args...>::type;

template <class Default, template <class...> class Op, class... Args>
using detected_or = traits_detail::detector<Default, void, Op, Args...>;

template <class Expected, template <typename...> class Op, class... Args>
using is_detected_exact = std::is_same<Expected, detected_t<Op, Args...>>;

////////////////////////////////////////////////////////////////////////////////
//                                                          [traits.conjunction]
//                                                          [traits.disjunction]
//                                                          [traits.negation]
#if defined(__cpp_lib_logical_traits)
using std::conjunction;
using std::disjunction;
using std::negation;
#elif defined(__cpp_lib_experimental_logical_traits)
using std::experimental::conjunction;
using std::experimental::disjunction;
using std::experimental::negation;
#else
template <class... Bools> using conjunction = meta::and_c<Bools::value...>;
template <class... Bools> using disjunction = meta::or_c<Bools::value...>;
template <class Bool> using negation = meta::not_c<Bool::value>;
#endif

////////////////////////////////////////////////////////////////////////////////
//                                                               [meta.constant]

using index = std::ptrdiff_t;
static constexpr index dyn{-1};

template <bool Bool> using bool_constant = std::integral_constant<bool, Bool>;
template <index I> using index_constant = std::integral_constant<index, I>;
template <int I> using int_constant = std::integral_constant<int, I>;
template <int I> static constexpr int_constant<I> int_c{};
template <bool B> static constexpr bool_constant<B> bool_c{};
template <index I> static constexpr index_constant<I> index_c{};

template <index N>
using make_index_sequence = std::make_integer_sequence<index, N>;

template <int N> using make_int_sequence = std::make_integer_sequence<int, N>;

////////////////////////////////////////////////////////////////////////////////
//                                                         [traits.is_invocable]
//                                                        [traits.invoke_result]

#ifdef __cpp_lib_is_invocable
using std::invoke_result;
using std::invoke_result_t;
using std::is_invocable;
using std::is_invocable_v;
using std::is_nothrow_invocable;
using std::is_nothrow_invocable_v;
#else
template <typename F, typename... Args>
using invoke_result_t =
    decltype(ranges::invoke(std::declval<F>(), std::declval<Args>()...));

template <typename F, typename... Args>
using is_invocable = is_detected<invoke_result_t, F, Args...>;
#endif

////////////////////////////////////////////////////////////////////////////////
//                                                                   [meta.size]

template <typename> struct list_size;
template <template <typename...> class L, typename... Vars>
struct list_size<L<Vars...>>
    : std::integral_constant<std::size_t, sizeof...(Vars)> {};
template <typename V>
static constexpr std::size_t list_size_v = list_size<V>::value;

////////////////////////////////////////////////////////////////////////////////
//                                                                 [meta.concat]

template <typename...> struct concat;
template <typename... Ls> using concat_t = typename concat<Ls...>::type;
template <typename L, typename... Ls> struct concat<L, Ls...> {
  using type = concat_t<L, concat_t<Ls...>>;
};
template <template <typename...> class L, typename... As, typename... Bs>
struct concat<L<As...>, L<Bs...>> {
  using type = L<As..., Bs...>;
};

////////////////////////////////////////////////////////////////////////////////
//                                                                   [meta.tail]

/// @brief Transforms any template list of type L<V1, V2, ...> to L<V2, ...>.
template <typename> struct tail;
template <template <typename...> class L, typename V, typename... Vars>
struct tail<L<V, Vars...>> {
  using type = L<Vars...>;
};
template <typename V> using tail_t = typename tail<V>::type;

////////////////////////////////////////////////////////////////////////////////
//                                                                   [meta.head]

/// @brief Returns the first type of a type list: L<V1, V2, ...> -> V1.
template <typename> struct head;
template <template <typename...> class L, typename V, typename... Vars>
struct head<L<V, Vars...>> {
  using type = V;
};
template <typename V> using head_t = typename head<V>::type;

////////////////////////////////////////////////////////////////////////////////
//                                                              [meta.list_cast]
//                                                               [meta.as_tuple]
//                                                          [meta.as_quantities]

/// @brief Transforms a type list `From` into a type list `To`.
template <template <typename...> class, typename> struct list_cast;
template <template <typename...> class To, template <typename...> class From,
          typename... Ts>
struct list_cast<To, From<Ts...>> {
  using type = To<Ts...>;
};
template <template <typename...> class To, typename From>
using list_cast_t = typename list_cast<To, From>::type;

/// @brief Transforms a type list to a std::tuple  L<Ts...> -> std::tuple<Ts...>
template <typename List> struct as_tuple : list_cast<std::tuple, List> {};
template <typename T, T... Is>
struct as_tuple<std::integer_sequence<T, Is...>> {
  using type = std::tuple<std::integral_constant<T, Is>...>;
};
template <typename List> using as_tuple_t = typename as_tuple<List>::type;

////////////////////////////////////////////////////////////////////////////////
//                                                             [meta.value_type]

template <typename Var> struct value_type {
  using type = typename Var::value_type;
};
template <typename Var> using value_type_t = typename value_type<Var>::type;

using ranges::data;
using ranges::size;

} // namespace fub

#endif // !TYPE_TRAITS_HPP
