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

#ifndef FUB_GRID_PATCH_VIEW_HPP
#define FUB_GRID_PATCH_VIEW_HPP

#include "fub/algorithm.hpp"
#include "fub/functional.hpp"
#include "fub/mdspan.hpp"
#include "fub/patch.hpp"
#include "fub/span.hpp"
#include "fub/tuple.hpp"
#include "fub/variables.hpp"

#include <range/v3/iterator_range.hpp>
#include <range/v3/utility/basic_iterator.hpp>
#include <range/v3/view/zip.hpp>

#include <unsupported/Eigen/CXX11/Tensor>

namespace fub {
template <typename, typename...> class patch_view;

namespace detail {
template <typename E> struct make_row_extents_fn;

template <index E, index... Es> struct make_row_extents_fn<extents<E, Es...>> {
  constexpr extents<E> operator()(index) const noexcept { return extents<E>(); }
};

template <index... Es> struct make_row_extents_fn<extents<dyn, Es...>> {
  constexpr extents<dyn> operator()(index len) const noexcept {
    return extents<dyn>(len);
  }
};

template <typename E>
constexpr auto make_row_extents(const E&, index len) noexcept {
  return make_row_extents_fn<E>()(len);
}

template <typename E> using row_extents_t = decltype(make_row_extents(E{}, 0));

template <index... Values>
constexpr std::ptrdiff_t covolume(const extents<Values...>& e) noexcept {
  constexpr int rank = extents<Values...>::rank;
  auto array = as_array(e);
  return accumulate(array.begin() + 1, array.end(), index(1),
                    std::multiplies<>());
}

template <typename E, typename... Vars> class row_cursor {
  std::tuple<typename variable_traits<Vars>::pointer...> m_pointers{};
  std::ptrdiff_t m_row_length{0};
  std::ptrdiff_t m_position{0};

public:
  struct mixin : ranges::basic_mixin<row_cursor> {
    mixin() = default;
    using ranges::basic_mixin<row_cursor>::basic_mixin;
    mixin(const patch_view<E, Vars...>& parent, std::ptrdiff_t pos = 0)
        : mixin(row_cursor(parent, pos)) {}
  };

  row_cursor() = default;

  explicit row_cursor(const patch_view<E, Vars...>& parent,
                      std::ptrdiff_t pos = 0)
      : m_pointers{parent.m_pointers}, m_row_length{parent.extents().get(0)},
        m_position{pos} {
    std::ptrdiff_t n = m_row_length * pos;
    for_each_tuple_element([n](auto& pointers) { std::advance(pointers, n); },
                           m_pointers);
  }

  patch_view<row_extents_t<E>, Vars...> read() const noexcept {
    const auto e = make_row_extents(E{}, m_row_length);
    const auto args = std::tuple_cat(std::make_tuple(e), m_pointers);
    return fub::make_from_tuple<patch_view<row_extents_t<E>, Vars...>>(args);
  }

  void next() noexcept {
    fub::for_each_tuple_element(
        [row_length = m_row_length](auto& pointers) {
          std::advance(pointers, row_length);
        },
        m_pointers);
    ++m_position;
  }

  bool equal(const row_cursor& other) const noexcept {
    return m_pointers == other.m_pointers && m_position == other.m_position;
  }
};

template <typename E, typename... Vars>
using row_iterator = ranges::basic_iterator<row_cursor<E, Vars...>>;

template <typename... Vars> class view_cursor {
  std::tuple<typename variable_traits<Vars>::pointer...> m_pointers{};
  std::ptrdiff_t m_position{0};

public:
  struct mixin : ranges::basic_mixin<view_cursor> {
    mixin() = default;
    using ranges::basic_mixin<view_cursor>::basic_mixin;
    template <typename E>
    mixin(const patch_view<E, Vars...>& parent, std::ptrdiff_t pos = 0)
        : mixin(view_cursor(parent, pos)) {}
  };

  view_cursor() = default;

  template <typename E>
  explicit view_cursor(const patch_view<E, Vars...>& parent,
                       std::ptrdiff_t pos = 0)
      : m_pointers{parent.m_pointers}, m_position{pos} {
    fub::for_each_tuple_element(
        [n = pos](auto& pointer) { std::advance(pointer, n); }, m_pointers);
  }

  quantities_ref<Vars...> read() const noexcept {
    return fub::make_from_tuple<quantities_ref<Vars...>>(m_pointers);
  }

  void next() noexcept {
    fub::for_each_tuple_element([](auto& pointer) { std::advance(pointer, 1); },
                                m_pointers);
    ++m_position;
  }

  bool equal(const view_cursor& other) const noexcept {
    return m_pointers == other.m_pointers && m_position == other.m_position;
  }
};

template <typename... Vars>
using patch_view_iterator = ranges::basic_iterator<view_cursor<Vars...>>;

template <typename Extents> class base_patch_view : public Extents {
public:
  base_patch_view() = default;
  base_patch_view(const Extents& e) : Extents(e) {}

  const Extents& extents() const noexcept { return *this; }
};

} // namespace detail

/// \brief span over multiple variables with uniform extents and layout mapping.
template <typename Extents, typename... Vars>
class patch_view : public detail::base_patch_view<Extents> {
  using base_type = detail::base_patch_view<Extents>;
  std::tuple<typename variable_traits<Vars>::pointer...> m_pointers{};

  friend class detail::row_cursor<Extents, Vars...>;
  friend class detail::view_cursor<Vars...>;

public:
  using extents_type = Extents;
  using variables = std::tuple<Vars...>;

  patch_view() = default;

  patch_view(const Extents& e,
             typename variable_traits<Vars>::pointer... pointers)
      : base_type(e), m_pointers{pointers...} {}

  template <typename Patch,
            typename = std::enable_if_t<is_patch<std::decay_t<Patch>>::value>>
  patch_view(Patch&& block)
      : base_type(block.extents()),
        m_pointers{block.template get<std::remove_const_t<Vars>>().data()...} {}

  template <typename Var> auto get() const noexcept {
    using Index = variable_find_index<std::remove_cv_t<Var>,
                                      meta::list<std::remove_cv_t<Vars>...>>;
    static_assert(Index::value < sizeof...(Vars),
                  "Could not find requested variable.");
    using ValueT = typename variable_traits<std::decay_t<Var>>::value_type;
    using VarInPack = std::tuple_element_t<Index::value, variables>;
    using ElementT = std::conditional_t<std::is_const<VarInPack>::value,
                                        const ValueT, ValueT>;
    return mdspan<ElementT, Extents>{std::get<Index::value>(m_pointers),
                                     this->extents()};
  }

  std::ptrdiff_t size() const noexcept { return this->extents().size(); }

  template <typename Var> auto operator[](Var) const noexcept {
    return get<Var>();
  }

  template <std::size_t N>
  quantities_ref<Vars...> operator()(const array<std::ptrdiff_t, N>& idx) const
      noexcept {
    return quantities_ref<Vars...>{
        std::addressof(fub::apply(get<Vars>(), idx))...};
  }

  template <typename... IndexTs>
  quantities_ref<Vars...> operator()(IndexTs... is) const noexcept {
    return this->operator()(array<std::ptrdiff_t, sizeof...(IndexTs)>{{is...}});
  }

  auto rows() const noexcept {
    detail::row_iterator<Extents, Vars...> first{*this};
    detail::row_iterator<Extents, Vars...> last{*this,
                                                detail::covolume(this->extents())};
    return ranges::make_iterator_range(first, last);
  }
};

/// \brief span over multiple variables with uniform extents and layout mapping.
template <int Size, typename... Vars>
class patch_view<extents<Size>, Vars...> : public detail::base_patch_view<extents<Size>> {
  using base_type = detail::base_patch_view<extents<Size>>;
  std::tuple<typename variable_traits<Vars>::pointer...> m_pointers{};

  friend class detail::row_cursor<extents<Size>, Vars...>;
  friend class detail::view_cursor<Vars...>;

public:
  using extents_type = ::fub::extents<Size>;
  using variables = std::tuple<Vars...>;

  patch_view() = default;

  patch_view(const extents_type& e,
             typename variable_traits<Vars>::pointer... pointers)
      : base_type(e), m_pointers{pointers...} {}

  template <typename Patch,
            typename = std::enable_if_t<is_patch<std::decay_t<Patch>>::value>>
  patch_view(Patch&& block)
      : base_type(block.extents()),
        m_pointers{block.template get<std::remove_const_t<Vars>>().data()...} {}

  template <typename Var> auto get() const noexcept {
    using Index = variable_find_index<std::remove_cv_t<Var>,
                                      meta::list<std::remove_cv_t<Vars>...>>;
    static_assert(Index::value < sizeof...(Vars),
                  "Could not find requested variable.");
    using ValueT = typename variable_traits<std::decay_t<Var>>::value_type;
    using VarInPack = std::tuple_element_t<Index::value, variables>;
    using ElementT = std::conditional_t<std::is_const<VarInPack>::value,
                                        const ValueT, ValueT>;
    return span<ElementT, Size>{std::get<Index::value>(m_pointers), Size};
  }

  std::ptrdiff_t size() const noexcept { return this->extents().size(); }

  template <typename Var> auto operator[](Var) const noexcept {
    return get<Var>();
  }

  template <std::size_t N>
  quantities_ref<Vars...> operator()(const array<std::ptrdiff_t, N>& idx) const
      noexcept {
    return quantities_ref<Vars...>{
        std::addressof(fub::apply(get<Vars>(), idx))...};
  }

  template <typename... IndexTs>
  quantities_ref<Vars...> operator()(IndexTs... is) const noexcept {
    return this->operator()(array<std::ptrdiff_t, sizeof...(IndexTs)>{{is...}});
  }

  quantities_ref<Vars...> first() const noexcept { return this->operator()(0); }

  quantities_ref<Vars...> last() const noexcept {
    return this->operator()(size() - 1);
  }

  auto begin() const noexcept {
    return detail::patch_view_iterator<Vars...>{*this, 0};
  }

  auto end() const noexcept {
    return detail::patch_view_iterator<Vars...>{*this, size()};
  }

  auto rows() const noexcept {
    detail::row_iterator<extents_type, Vars...> first{*this};
    detail::row_iterator<extents_type, Vars...> last{
        *this, detail::covolume(this->extents())};
    return ranges::make_iterator_range(first, last);
  }
};

template <typename Extents, typename... Vars>
using const_patch_view = patch_view<Extents, std::add_const_t<Vars>...>;

template <typename T> struct view_add_const;
template <typename T> using view_add_const_t = typename view_add_const<T>::type;
template <typename Extents, typename... Vars>
struct view_add_const<patch_view<Extents, Vars...>> {
  using type = patch_view<Extents, std::add_const_t<Vars>...>;
};

template <typename Extents, typename Storage, template <typename...> class V,
          typename... Vars>
auto make_view(patch<V<Vars...>, Extents, Storage>& p) noexcept {
  return patch_view<Extents, Vars...>(p);
}

template <typename Extents, typename Storage, template <typename...> class V,
          typename... Vars>
auto make_view(const patch<V<Vars...>, Extents, Storage>& p) noexcept {
  return patch_view<Extents, std::add_const_t<Vars>...>(p);
}

template <typename Vars, index Size> struct view_row;
template <typename V, index Size>
using view_row_t = typename view_row<V, Size>::type;
template <template <typename...> class V, index Size, typename... Vs>
struct view_row<V<Vs...>, Size> {
  using type = patch_view<extents<Size>, Vs...>;
};

template <typename E, typename... Vars>
struct add_flux<patch_view<E, Vars...>> {
  using type = patch_view<E, add_flux_t<Vars>...>;
};

template <typename E, typename... Vars>
struct remove_flux<patch_view<E, Vars...>> {
  using type = patch_view<E, remove_flux_t<Vars>...>;
};

template <typename patch_view> struct is_view : std::false_type {};
template <typename T> static constexpr bool is_view_v = is_view<T>::value;

template <typename E, typename... Vars>
struct is_view<patch_view<E, Vars...>> : std::true_type {};

template <typename T, int Dim> struct view_static_extent;
template <int Dim, typename E, typename... Vs>
struct view_static_extent<patch_view<E, Vs...>, Dim>
    : std::integral_constant<std::ptrdiff_t, E().get(Dim)> {};

template <typename T, int Dim>
static constexpr std::ptrdiff_t view_static_extent_v =
    view_static_extent<T, Dim>::value;

template <typename V> struct view_rank;
template <typename E, typename... Vs>
struct view_rank<patch_view<E, Vs...>> : std::integral_constant<int, E::rank> {
};
template <typename V> static constexpr int view_rank_v = view_rank<V>::value;

template <typename V> struct view_variables;
template <typename E, typename... Vs>
struct view_variables<patch_view<E, Vs...>> {
  using type = std::tuple<Vs...>;
};
template <typename V> using view_variables_t = typename view_variables<V>::type;

template <typename F, typename... Views>
std::enable_if_t<conjunction<is_view<Views>...>::value, F>
for_each_row(F function, const Views&... views) {
  auto all_rows = ranges::view::zip(views.rows()...);
  for (auto&& rows : all_rows) {
    fub::apply(std::ref(function), rows);
  }
  return function;
}

////////////////////////////////////////////////////////////////////////////////
// drop / take with static extents

template <index N, index Len, typename... Vars>
auto drop(const patch_view<extents<Len>, Vars...>& view) noexcept {
  return patch_view<extents<Len - N>, Vars...>(
      extents<Len - N>(), drop<N>(view.template get<Vars>()).data()...);
}

template <index N, index Len, typename... Vars>
auto take(const patch_view<extents<Len>, Vars...>& view) noexcept {
  return patch_view<extents<N>, Vars...>(
      extents<N>(), take<N>(view.template get<Vars>()).data()...);
}

template <index N, index Len, typename... Vars>
auto rtake(const patch_view<extents<Len>, Vars...>& view) noexcept {
  return drop<Len - N>(view);
}

template <index N, index Len, typename... Vars>
auto rdrop(const patch_view<extents<Len>, Vars...>& view) noexcept {
  return take<Len - N>(view);
}

template <typename View> struct view_variable_tuple;
template <typename E, typename... Vars>
struct view_variable_tuple<patch_view<E, Vars...>> {
  using type = std::tuple<Vars...>;
};
template <typename View>
using view_variable_tuple_t = typename view_variable_tuple<View>::type;

template <typename View> struct view_extents;
template <typename E, typename... Vars>
struct view_extents<patch_view<E, Vars...>> {
  using type = E;
};
template <typename View>
using view_extents_t = typename view_extents<View>::type;

template <typename... Vars> struct view_join_variables;

template <typename E, typename... Vars, typename... Ts>
struct view_join_variables<patch_view<E, Vars...>, Ts...> {
  using type = std::tuple<std::remove_const_t<Vars>...>;
};

template <typename... Vars>
using view_join_variables_t = typename view_join_variables<Vars...>::type;

template <typename... RowViews> auto join(const RowViews&... rows) {
  static_assert(conjunction<is_view<RowViews>...>::value,
                "Only views can be joined.");
  static_assert(
      conjunction<
          bool_constant<view_extents_t<RowViews>::rank_dynamic == 0>...>::value,
      "Only statically sized views can be joined.");
  static_assert(
      conjunction<bool_constant<view_extents_t<RowViews>::rank == 1>...>::value,
      "Only one-dimensional views an be joined.");
  static constexpr array<index, sizeof...(RowViews)> sizes{
      {view_static_extent<RowViews, 0>::value...}};
  static constexpr int total_size = accumulate(sizes, index(0));
  patch<view_join_variables_t<RowViews...>, extents<total_size>,
        automatic_storage_descriptor>
      joined_row;
  auto view = make_view(joined_row);
  auto it = view.begin();
  (void)std::initializer_list<int>{
      ((void)(it = std::copy(rows.begin(), rows.end(), it)), 42)...};
  return joined_row;
}

/// \brief Loads all quantities which are refered by a specified view into simd
/// registers.
///
/// Use this to transform a view and position index into a simd vector when
/// doing vecotrised computations.
template <typename Abi, typename E, typename... Vs>
quantities<add_simd_t<std::remove_const_t<Vs>, Abi>...>
load(const patch_view<E, Vs...>& view, index shift = 0) noexcept {
  auto load = [=](auto x) {
    using T = typename variable_traits<std::decay_t<decltype(x)>>::value_type;
    simd<T, Abi> pack;
    pack.copy_from(x.data() + shift, element_alignment);
    return pack;
  };
  return {load(view.template get<Vs>())...};
}

template <typename E, typename... Vs>
void store(const quantities<Vs...>& q, const patch_view<E, Vs...>& view,
           index shift) noexcept {
  view(shift) = q;
}

template <typename Abi, typename E, typename... Vs>
void store(const add_simd_t<quantities<Vs...>, Abi>& q,
           const patch_view<E, Vs...>& view, index shift) noexcept {
  auto store = [=](const simd<double, Abi>& pack, double* pointer) {
    pack.copy_to(pointer + shift, element_alignment);
  };
  fub::for_each_tuple_element(
      [&](auto var) { store(q[var], view[var].data()); }, std::tuple<Vs...>{});
}

template <typename V, typename Abi> struct simd_proxy {
  using Variables = view_variables_t<V>;
  using simd_type = add_simd_t<as_quantities_t<Variables>, Abi>;

  V view;
  index shift;

  simd_proxy(V v, index s) noexcept : view{v}, shift{s} {}
  simd_proxy(const simd_proxy&);
  simd_proxy& operator=(const simd_proxy&);

  operator simd_type() const noexcept { return load<Abi>(view, shift); }

  void operator=(const simd_type& quantities) noexcept {
    store<Abi>(quantities, view, shift);
  }
};

template <typename Abi, typename E, typename... Vs>
std::enable_if_t<conjunction<std::is_const<Vs>...>::value,
                 quantities<add_simd_t<std::remove_const_t<Vs>, Abi>...>>
load_or_proxy(const patch_view<E, Vs...>& view, index shift = 0) noexcept {
  return load<Abi>(view, shift);
}

template <typename Abi, typename E, typename... Vs>
std::enable_if_t<conjunction<negation<std::is_const<Vs>>...>::value,
                 simd_proxy<patch_view<E, Vs...>, Abi>>
load_or_proxy(const patch_view<E, Vs...>& view, index shift = 0) noexcept {
  return simd_proxy<patch_view<E, Vs...>, Abi>{view, shift};
}

namespace detail {
template <typename Abi, typename F, typename FirstRow, typename... Rows>
std::enable_if_t<is_simd_abi_v<Abi>, F>
for_each_simd_impl(const Abi& abi, F f, const FirstRow& first,
                   const Rows&... rows);

template <typename F, typename Row, typename... Rows>
F for_each_simd_impl(index_constant<0>, F&& f, const Row&, const Rows&...) {
  return std::forward<F>(f);
}

template <index Size, typename F, typename Row, typename... Rows>
F for_each_simd_impl(index_constant<Size>, F f, const Row& row,
                     const Rows&... rows) {
  using Extents = view_extents_t<Row>;
  static constexpr index size = Extents().size();
  return for_each_simd_impl(simd_abi::scalar(), std::move(f),
                            drop<size - Size>(row), drop<size - Size>(rows)...);
}

template <typename Abi, typename F, typename FirstRow, typename... Rows>
std::enable_if_t<is_simd_abi_v<Abi>, F>
for_each_simd_impl(const Abi& abi, F f, const FirstRow& first,
                   const Rows&... rows) {
  static_assert(conjunction<is_view<FirstRow>, is_view<Rows>...>::value,
                "All parameters need to be views.");
  static_assert(conjunction<std::is_same<view_extents_t<FirstRow>,
                                         view_extents_t<Rows>>...>::value,
                "All parameters need to have the same size.");
  using Extents = view_extents_t<FirstRow>;
  static_assert(Extents::rank_dynamic == 0,
                "patch_view size needs to be a compile-time known.");
  static_assert(Extents::rank == 1,
                "Only one-dimensional views are supported.");
  static constexpr index width = simd<double, Abi>::size();
  static constexpr index size = Extents().size();
  static constexpr index simd_size = size / width;
  for (int simd_index = 0; simd_index < simd_size; ++simd_index) {
    std::ptrdiff_t shift = simd_index * width;
    fub::invoke(std::ref(f), abi, load_or_proxy<Abi>(first, shift),
                load_or_proxy<Abi>(rows, shift)...);
  }
  static constexpr index rest = size % width;
  return for_each_simd_impl(index_c<rest>, std::move(f), first, rows...);
}

} // namespace detail

template <typename Abi, typename F, typename FirstRow, typename... Rows>
std::enable_if_t<is_simd_abi_v<Abi>, F>
for_each_simd(const Abi& abi, F f, const FirstRow& first, const Rows&... rows) {
  return detail::for_each_simd_impl(abi, std::move(f), first, rows...);
}

template <typename F, typename... Rows>
std::enable_if_t<conjunction<is_view<Rows>...>::value, F>
for_each_simd(F f, const Rows&... rows) {
  return for_each_simd(simd_abi::native<double>(), std::move(f), rows...);
}

//                                                     [function.permutate]

template <int A, int B, index... Es>
constexpr array<index, sizeof...(Es)>
make_permutated_extent_array(const extents<Es...>&) {
  array<index, sizeof...(Es)> e{{Es...}};
  index eA = e[A];
  e[A] = e[B];
  e[B] = eA;
  return e;
}

template <int A, int B, typename E, std::size_t... Is>
constexpr auto permutate_extents_impl(const E&, std::index_sequence<Is...>) {
  constexpr array<index, sizeof...(Is)> es =
      make_permutated_extent_array<A, B>(E{});
  return extents<es[Is]...>();
}

template <int A, int B, index... Es>
constexpr auto permutate_extents(const extents<Es...>& e) {
  return permutate_extents_impl<A, B>(
      e, std::make_index_sequence<sizeof...(Es)>());
}

template <typename E, int A, int B>
using permutate_extents_t = decltype(permutate_extents<A, B>(E()));

template <typename S> struct tensor_map;
template <typename S> using tensor_map_t = typename tensor_map<S>::type;
template <typename T, index... Es>
struct tensor_map<mdspan<T, extents<Es...>>> {
  using type = Eigen::TensorMap<Eigen::Tensor<T, sizeof...(Es)>>;
};
template <typename T, index... Es>
tensor_map_t<mdspan<T, extents<Es...>>>
make_tensor(const mdspan<T, extents<Es...>>& view) {
  return tensor_map_t<mdspan<T, extents<Es...>>>(view.data(), Es...);
}

template <int A, int B, int Rank>
constexpr array<int, Rank> make_permutation_index() noexcept {
  auto index = fub::apply(
      [](auto... is) { return array<int, sizeof...(is)>{{is()...}}; },
      as_tuple_t<std::make_integer_sequence<int, Rank>>());
  std::swap(index[A], index[B]);
  return index;
}

template <int A, int B, typename V> struct variable_permutate {
  using type = V;
};
template <int A, int B, typename V>
using variable_permutate_t = typename variable_permutate<A, B, V>::type;
template <int A, int B, template <int> class V>
struct variable_permutate<A, B, V<A>> {
  using type = V<B>;
};
template <int A, int B, template <int> class V>
struct variable_permutate<A, B, V<B>> {
  using type = V<A>;
};

/// @brief Returns a patch which is a copy of the patch view but has permutated
/// axis.
///
/// This deeply copies all data and rearranged them in memory such that the data
/// in the specified axis will be accessible well.
template <int A, int B, typename ViewA, typename ViewB>
void permutate(const ViewA& dest, const ViewB& source) {
  for_each_tuple_element(
      [&](auto variable) {
        auto permutated = variable_permutate_t<A, B, decltype(variable)>();
        auto st = make_tensor(source[variable]);
        auto dt = make_tensor(dest[permutated]);
        auto index = make_permutation_index<A, B, view_rank_v<ViewA>>();
        dt = st.shuffle(index);
      },
      view_variables_t<ViewA>());
}

template <int A, int B, typename Extents, typename... Vars>
patch<std::tuple<std::remove_const_t<Vars>...>,
      permutate_extents_t<Extents, A, B>>
permutate(const patch_view<Extents, Vars...>& source) {
  patch<std::tuple<std::remove_const_t<Vars>...>,
        permutate_extents_t<Extents, A, B>>
      dest;
  permutate<A, B>(make_view(dest), source);
  return dest;
}

template <int Dim, int Width, typename ViewA, typename ViewB>
void slice_left(const ViewA& dest, const ViewB& source) {
  for_each_tuple_element(
      [&](auto variable) {
        auto st = make_tensor(source[variable]);
        auto dt = make_tensor(dest[variable]);
        array<index, view_rank_v<ViewA>> origin{};
        origin[Dim] = source.extents().get(Dim) - Width;
        array<index, view_rank_v<ViewA>> extent = dest.extents();
        dt = st.slice(origin, extent);
      },
      view_variables_t<ViewA>());
}

template <int Dim, int Width, typename Extents, typename... Vars>
auto slice_left(const patch_view<Extents, Vars...>& source) {
  auto changed_extents = replace_extent<Dim, Width>(source.extents());
  auto dest = make_patch(std::tuple<Vars...>(), changed_extents);
  slice_left<Dim, Width>(make_view(dest), source);
  return dest;
}

template <int Dim, int Width, typename ViewA, typename ViewB>
void slice_right(const ViewA& dest, const ViewB& source) {
  for_each_tuple_element(
      [&](auto variable) {
        auto st = make_tensor(source[variable]);
        auto dt = make_tensor(dest[variable]);
        array<index, view_rank_v<ViewA>> origin{};
        array<index, view_rank_v<ViewA>> extent = dest.extents();
        dt = st.slice(origin, extent);
      },
      view_variables_t<ViewA>());
}

template <int Dim, int Width, typename Extents, typename... Vars>
auto slice_right(const patch_view<Extents, Vars...>& source) {
  auto changed_extents = replace_extent<Dim, Width>(source.extents());
  auto dest = make_patch(std::tuple<Vars...>(), changed_extents);
  slice_right<Dim, Width>(make_view(dest), source);
  return dest;
}

} // namespace fub

#endif // !BLOCKVIEW_HPP
