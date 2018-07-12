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

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#endif
#include <unsupported/Eigen/CXX11/Tensor>
#ifdef __clang__
#pragma clang diagnostic pop
#endif

namespace fub {
template <typename, typename...> class patch_view;
template <index, typename...> class row_view;

namespace detail {
template <typename E> struct make_row_extents_fn;

template <index E, index... Es> struct make_row_extents_fn<extents<E, Es...>> {
  constexpr extents<E> operator()(index) const noexcept { return extents<E>(); }
};

template <index... Es>
struct make_row_extents_fn<extents<dynamic_extent, Es...>> {
  constexpr extents<dynamic_extent> operator()(index len) const noexcept {
    return extents<dynamic_extent>(len);
  }
};

template <typename E>
constexpr auto make_row_extents(const E&, index len) noexcept {
  return make_row_extents_fn<E>()(len);
}

template <typename E> struct row_extents;
template <index Size, index... Rest>
struct row_extents<extents<Size, Rest...>> : index_constant<Size> {};

template <index... Values>
constexpr std::ptrdiff_t covolume(const extents<Values...>& e) noexcept {
  auto array = as_array(e);
  return fub::accumulate(array.begin() + 1, array.end(), index(1),
                         std::multiplies<>());
}

template <typename E, typename... Vars> class row_cursor {
  template <typename V>
  using pointer = std::conditional_t<
      std::is_const<V>::value,
      typename variable_traits<std::remove_const_t<V>>::const_pointer,
      typename variable_traits<V>::pointer>;
  std::tuple<pointer<Vars>...> m_pointers{};
  std::ptrdiff_t m_row_length{0};
  std::ptrdiff_t m_position{0};

public:
  struct mixin : ranges::basic_mixin<row_cursor> {
    mixin() = default;
    using ranges::basic_mixin<row_cursor>::basic_mixin;
    mixin(const patch_view<E, Vars...>& parent, std::ptrdiff_t pos = 0)
        : mixin(row_cursor(parent, pos)) {}
    template <index N>
    mixin(const row_view<N, Vars...>& parent, std::ptrdiff_t pos = 0)
        : mixin(row_cursor(parent, pos)) {}
  };

  row_cursor() = default;

  explicit row_cursor(const patch_view<E, Vars...>& parent,
                      std::ptrdiff_t pos = 0)
      : m_pointers{parent.m_pointers}, m_row_length{parent.extents().get(0)},
        m_position{pos} {
    std::ptrdiff_t n = m_row_length * pos;
    for_each_tuple_element([n](auto& ps) { ps = ps + n; }, m_pointers);
  }

  template <index N>
  explicit row_cursor(const row_view<N, Vars...>& parent,
                      std::ptrdiff_t pos = 0)
      : m_pointers{parent.m_pointers}, m_row_length{parent.extents().size()},
        m_position{pos} {
    std::ptrdiff_t n = m_row_length * pos;
    for_each_tuple_element([n](auto& ps) { ps = ps + n; }, m_pointers);
  }

  row_view<row_extents<E>::value, Vars...> read() const noexcept {
    const auto e = make_row_extents(E{}, m_row_length);
    const auto args = std::tuple_cat(std::make_tuple(e), m_pointers);
    return fub::make_from_tuple<row_view<row_extents<E>::value, Vars...>>(args);
  }

  void next() noexcept {
    for_each_tuple_element(
        [row_length = m_row_length](auto& ps) { ps = ps + row_length; },
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
  template <typename V>
  using pointer = std::conditional_t<
      std::is_const<V>::value,
      typename variable_traits<std::remove_const_t<V>>::const_pointer,
      typename variable_traits<V>::pointer>;
  std::tuple<pointer<Vars>...> m_pointers{};
  std::ptrdiff_t m_position{0};

public:
  struct mixin : ranges::basic_mixin<view_cursor> {
    mixin() = default;
    using ranges::basic_mixin<view_cursor>::basic_mixin;
    template <typename E>
    mixin(const patch_view<E, Vars...>& parent, std::ptrdiff_t pos = 0)
        : mixin(view_cursor(parent, pos)) {}

    template <index E>
    mixin(const row_view<E, Vars...>& parent, std::ptrdiff_t pos = 0)
        : mixin(view_cursor(parent, pos)) {}
  };

  view_cursor() = default;

  template <typename E>
  explicit view_cursor(const patch_view<E, Vars...>& parent,
                       std::ptrdiff_t pos = 0)
      : m_pointers{parent.m_pointers}, m_position{pos} {
    fub::for_each_tuple_element([n = pos](auto& ptr) { ptr = ptr + n; },
                                m_pointers);
  }

  template <index E>
  explicit view_cursor(const row_view<E, Vars...>& parent,
                       std::ptrdiff_t pos = 0)
      : m_pointers{parent.m_pointers}, m_position{pos} {
    fub::for_each_tuple_element([n = pos](auto& ptr) { ptr = ptr + n; },
                                m_pointers);
  }

  quantities_ref<Vars...> read() const noexcept {
    return fub::make_from_tuple<quantities_ref<Vars...>>(m_pointers);
  }

  void next() noexcept {
    fub::for_each_tuple_element([](auto& ptr) { ptr = ptr + index(1); },
                                m_pointers);
    ++m_position;
  }

  bool equal(const view_cursor& other) const noexcept {
    return m_pointers == other.m_pointers && m_position == other.m_position;
  }
};

template <typename... Vars>
using patch_view_iterator = ranges::basic_iterator<view_cursor<Vars...>>;

template <typename Extents> class patch_view_base : public Extents {
public:
  patch_view_base() = default;
  patch_view_base(const Extents& e) : Extents(e) {}

  const Extents& extents() const noexcept { return *this; }
};

} // namespace detail

/// \brief span over multiple variables with uniform extents and layout mapping.
template <typename Extents, typename... Vars>
class patch_view : public detail::patch_view_base<Extents> {
  using base_type = detail::patch_view_base<Extents>;
  template <typename V>
  using pointer = std::conditional_t<
      std::is_const<V>::value,
      typename variable_traits<std::remove_const_t<V>>::const_pointer,
      typename variable_traits<V>::pointer>;
  std::tuple<pointer<Vars>...> m_pointers{};

  friend class detail::row_cursor<Extents, Vars...>;
  friend class detail::view_cursor<Vars...>;

public:
  using extents_type = Extents;
  using variables = std::tuple<Vars...>;

  patch_view() = default;

  patch_view(const Extents& e, pointer<Vars>... ptr)
      : base_type(e), m_pointers{ptr...} {}

  template <typename Patch,
            typename = std::enable_if_t<is_patch<std::decay_t<Patch>>::value>>
  patch_view(Patch&& block)
      : base_type(block.extents()),
        m_pointers{variable_traits<std::remove_const_t<Vars>>::get_pointer(
            block.template get<std::remove_const_t<Vars>>())...} {}

  template <typename Var> auto get() const noexcept {
    using Index = variable_find_index<Var, meta::list<Vars...>>;
    static_assert(Index::value < sizeof...(Vars),
                  "Could not find requested variable.");
    using Holder =
        std::remove_const_t<std::tuple_element_t<Index::value, variables>>;
    return variable_traits<Holder>::view(
        Var(), std::get<Index::value>(m_pointers), this->extents());
  }

  std::ptrdiff_t size() const noexcept { return this->extents().size(); }

  template <typename Var> auto operator[](Var) const noexcept {
    return get<Var>();
  }

  template <std::size_t N>
  quantities_ref<Vars...>
  operator()(const std::array<std::ptrdiff_t, N>& idx) const noexcept {
    layout_left::mapping<extents_type> mapping;
    auto offset = fub::apply(mapping, idx);
    return fub::apply(
        [=](auto... Is) {
          return quantities_ref<Vars...>{
              (std::get<decltype(Is)::value>(m_pointers) + offset)...};
        },
        as_tuple_t<decltype(make_index_sequence<sizeof...(Vars)>())>{});
  }

  template <typename... IndexTs>
  quantities_ref<Vars...> operator()(IndexTs... is) const noexcept {
    return this->operator()(
        std::array<std::ptrdiff_t, sizeof...(IndexTs)>{{is...}});
  }

  auto rows() const noexcept {
    detail::row_iterator<Extents, Vars...> first{*this};
    detail::row_iterator<Extents, Vars...> last{
        *this, detail::covolume(this->extents())};
    return ranges::make_iterator_range(first, last);
  }
};

/// \brief span over multiple variables with uniform extents and layout mapping.
template <index Size, typename... Vars>
class row_view : public detail::patch_view_base<extents<Size>> {
public:
  static_assert(Size == dynamic_extent || Size > 0, "Invalid Row Size");

  using base_type = detail::patch_view_base<extents<Size>>;

  template <typename V>
  using pointer = std::conditional_t<
      std::is_const<V>::value,
      typename variable_traits<std::remove_const_t<V>>::const_pointer,
      typename variable_traits<V>::pointer>;
  std::tuple<pointer<Vars>...> m_pointers;

  friend class detail::row_cursor<extents<Size>, Vars...>;
  friend class detail::view_cursor<Vars...>;

public:
  using extents_type = ::fub::extents<Size>;
  using variables = std::tuple<Vars...>;

  row_view() = default;

  row_view(const extents_type& e, pointer<Vars>... ps)
      : base_type(e), m_pointers{ps...} {}

  template <typename Patch,
            typename = std::enable_if_t<is_patch<std::decay_t<Patch>>::value>>
  row_view(Patch&& block)
      : base_type(block.extents()),
        m_pointers{variable_traits<std::remove_const_t<Vars>>::get_pointer(
            block.template get<std::remove_const_t<Vars>>())...} {}

  template <typename Var> auto get() const noexcept {
    using Index = variable_find_index<Var, meta::list<Vars...>>;
    static_assert(Index::value < sizeof...(Vars),
                  "Could not find requested variable.");
    using Holder =
        std::remove_const_t<std::tuple_element_t<Index::value, variables>>;
    auto mdspan = variable_traits<Holder>::view(
        Var(), std::get<Index::value>(m_pointers), this->extents());
    return mdspan.span();
  }

  index size() const noexcept { return this->extents().size(); }

  template <typename Var> auto operator[](Var) const noexcept {
    return get<Var>();
  }

  quantities_ref<Vars...>
  operator()(const std::array<std::ptrdiff_t, 1>& idx) const noexcept {
    index offset = idx[0];
    return fub::apply(
        [=](auto... Is) {
          return quantities_ref<Vars...>{
              (std::get<decltype(Is)::value>(m_pointers) + offset)...};
        },
        as_tuple_t<decltype(make_index_sequence<sizeof...(Vars)>())>{});
  }

  template <typename... IndexTs>
  quantities_ref<Vars...> operator()(IndexTs... is) const noexcept {
    return this->operator()(
        std::array<std::ptrdiff_t, sizeof...(IndexTs)>{{is...}});
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
    detail::row_iterator<extents<Size>, Vars...> first{*this};
    detail::row_iterator<extents<Size>, Vars...> last{
        *this, detail::covolume(this->extents())};
    return ranges::make_iterator_range(first, last);
  }
};

////////////////////////////////////////////////////////////////////////////////
// drop / take with static extents

template <index N, index Len, typename... Vs,
          typename = std::enable_if_t<(Len > 0) && (N < Len)>>
auto drop(const row_view<Len, Vs...>& view) noexcept {
  return fub::apply(
      [](auto... pointers) {
        return row_view<Len - N, Vs...>(extents<Len - N>(), (pointers + N)...);
      },
      view.m_pointers);
}

template <index N, typename... Vs>
auto drop(const row_view<dynamic_extent, Vs...>& view) noexcept {
  const index size = view.extents().size() - N;
  assert(size > 0);
  return fub::apply(
      [=](auto... pointers) {
        return row_view<dynamic_extent, Vs...>(extents<dynamic_extent>(size),
                                               (pointers + N)...);
      },
      view.m_pointers);
}

template <index N, index Len, typename... Vs,
          typename = std::enable_if_t<(Len > 0) && (N <= Len)>>
auto take(const row_view<Len, Vs...>& view) noexcept {
  return fub::apply(
      [](auto... pointers) {
        return row_view<N, Vs...>(extents<N>(), pointers...);
      },
      view.m_pointers);
}

template <index N, typename... Vs>
auto take(const row_view<dynamic_extent, Vs...>& view) noexcept {
  assert(N <= view.extents().size());
  return fub::apply(
      [](auto... pointers) {
        return row_view<N, Vs...>(extents<N>(), pointers...);
      },
      view.m_pointers);
}

template <index N, index Len, typename... Vs,
          typename = std::enable_if_t<(Len > 0) && (N <= Len)>>
auto rtake(const row_view<Len, Vs...>& view) noexcept {
  return drop<Len - N>(view);
}

template <index N, typename... Vs>
auto rtake(const row_view<dynamic_extent, Vs...>& view) noexcept {
  const index offset = view.extents().size() - N;
  assert(offset > 0);
  return fub::apply(
      [=](auto... pointers) {
        return row_view<N, Vs...>(extents<N>(), (pointers + offset)...);
      },
      view.m_pointers);
}

template <typename... Vs>
auto rtake(const row_view<dynamic_extent, Vs...>& view, index n) noexcept {
  const index offset = view.extents().size() - n;
  assert(offset > 0);
  return fub::apply(
      [=](auto... pointers) {
        return row_view<dynamic_extent, Vs...>(extents<dynamic_extent>(n),
                                               (pointers + offset)...);
      },
      view.m_pointers);
}

template <index N, index Len, typename... Vs,
          typename = std::enable_if_t<(Len > 0) && (N < Len)>>
auto rdrop(const row_view<Len, Vs...>& view) noexcept {
  return take<Len - N>(view);
}

template <index N, typename... Vs>
auto rdrop(const row_view<dynamic_extent, Vs...>& view) noexcept {
  const index size = view.extents().size() - N;
  assert(size > 0);
  return fub::apply(
      [=](auto... pointers) {
        return row_view<dynamic_extent, Vs...>(extents<dynamic_extent>(size),
                                               pointers...);
      },
      view.m_pointers);
}

template <typename Extents, typename... Vars>
using const_patch_view = patch_view<Extents, std::add_const_t<Vars>...>;

template <typename T> struct view_add_const;
template <typename T> using view_add_const_t = typename view_add_const<T>::type;
template <typename Extents, typename... Vars>
struct view_add_const<patch_view<Extents, Vars...>> {
  using type = patch_view<Extents, std::add_const_t<Vars>...>;
};
template <std::ptrdiff_t Extents, typename... Vars>
struct view_add_const<row_view<Extents, Vars...>> {
  using type = row_view<Extents, std::add_const_t<Vars>...>;
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

template <index E, typename Storage, template <typename...> class V,
          typename... Vars>
auto make_view(patch<V<Vars...>, extents<E>, Storage>& p) noexcept {
  return row_view<E, Vars...>(p);
}

template <index E, typename Storage, template <typename...> class V,
          typename... Vars>
auto make_view(const patch<V<Vars...>, extents<E>, Storage>& p) noexcept {
  return row_view<E, std::add_const_t<Vars>...>(p);
}

template <typename T>
using patch_view_t = decltype(make_view(std::declval<T>()));

template <typename Vars, index Size> struct view_row;
template <typename V, index Size>
using view_row_t = typename view_row<V, Size>::type;
template <template <typename...> class V, index Size, typename... Vs>
struct view_row<V<Vs...>, Size> {
  using type = row_view<Size, Vs...>;
};

template <typename E, typename... Vars>
struct add_flux<patch_view<E, Vars...>> {
  using type = patch_view<E, add_flux_t<Vars>...>;
};

template <typename E, typename... Vars>
struct remove_flux<patch_view<E, Vars...>> {
  using type = patch_view<E, remove_flux_t<Vars>...>;
};

template <typename T> struct is_view : std::false_type {};
template <typename T> static constexpr bool is_view_v = is_view<T>::value;

template <typename E, typename... Vars>
struct is_view<patch_view<E, Vars...>> : std::true_type {};
template <std::ptrdiff_t Size, typename... Vars>
struct is_view<row_view<Size, Vars...>> : std::true_type {};

template <typename T, int Dim> struct view_static_extent;
template <int Dim, typename E, typename... Vs>
struct view_static_extent<patch_view<E, Vs...>, Dim>
    : std::integral_constant<std::ptrdiff_t, E().get(Dim)> {};
template <std::ptrdiff_t Size, typename... Vs>
struct view_static_extent<row_view<Size, Vs...>, 0>
    : std::integral_constant<std::ptrdiff_t, Size> {};

template <typename T, int Dim>
static constexpr std::ptrdiff_t view_static_extent_v =
    view_static_extent<T, Dim>::value;

template <typename V> struct view_rank;
template <typename E, typename... Vs>
struct view_rank<patch_view<E, Vs...>> : std::integral_constant<int, E::rank> {
};
template <index E, typename... Vs>
struct view_rank<row_view<E, Vs...>> : std::integral_constant<int, 1> {};
template <typename V> static constexpr int view_rank_v = view_rank<V>::value;

template <typename V> struct view_variables;
template <typename E, typename... Vs>
struct view_variables<patch_view<E, Vs...>> {
  using type = flatten_variables_t<Vs...>;
};
template <index E, typename... Vs> struct view_variables<row_view<E, Vs...>> {
  using type = flatten_variables_t<Vs...>;
};
template <typename V> using view_variables_t = typename view_variables<V>::type;

template <typename F, typename... Views>
std::enable_if_t<conjunction<is_view<Views>...>::value, F>
for_each_row(F function, const Views&... views) {
  auto all_rows = ranges::view::zip(views.rows()...);
  for (auto&& rows : all_rows) {
    fub::apply(function, rows);
  }
  return function;
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
template <index E, typename... Vars> struct view_extents<row_view<E, Vars...>> {
  using type = extents<E>;
};
template <typename View>
using view_extents_t = typename view_extents<View>::type;

template <typename... Vars> struct view_join_variables;

template <typename E, typename... Vars, typename... Ts>
struct view_join_variables<patch_view<E, Vars...>, Ts...> {
  using type = std::tuple<std::remove_const_t<Vars>...>;
};
template <index E, typename... Vars, typename... Ts>
struct view_join_variables<row_view<E, Vars...>, Ts...> {
  using type = std::tuple<std::remove_const_t<Vars>...>;
};

template <typename... Vars>
using view_join_variables_t = typename view_join_variables<Vars...>::type;

namespace detail {
template <typename... RowViews>
auto join_impl(std::true_type, const RowViews&... rows) {
  static_assert(conjunction<is_view<RowViews>...>::value,
                "Only views can be joined.");
  static_assert(
      conjunction<
          bool_constant<view_extents_t<RowViews>::rank_dynamic() == 0>...>::value,
      "Only statically sized views can be joined.");
  static_assert(
      conjunction<bool_constant<view_extents_t<RowViews>::rank() == 1>...>::value,
      "Only one-dimensional views an be joined.");
  static constexpr index sizes[sizeof...(RowViews)]{
      view_static_extent<RowViews, 0>::value...};
  static constexpr index total_size = fub::accumulate(sizes, index(0));
  patch<view_join_variables_t<RowViews...>, extents<total_size>> joined_row;
  auto view = make_view(joined_row);
  auto it = view.begin();
  (void)std::initializer_list<int>{
      ((void)(it = std::copy(rows.begin(), rows.end(), it)), 42)...};
  return joined_row;
}
template <typename... RowViews>
auto join_impl(std::false_type, const RowViews&... rows) {
  static_assert(conjunction<is_view<RowViews>...>::value,
                "Only views can be joined.");
  static_assert(
      conjunction<bool_constant<view_extents_t<RowViews>::rank() == 1>...>::value,
      "Only one-dimensional views an be joined.");
  const std::array<index, sizeof...(RowViews)> sizes{
      {rows.extents().get(0)...}};
  const extents<dynamic_extent> total_size{fub::accumulate(sizes, index(0))};
  patch<view_join_variables_t<RowViews...>, extents<dynamic_extent>> joined_row(
      total_size);
  auto view = make_view(joined_row);
  auto it = view.begin();
  (void)std::initializer_list<int>{
      ((void)(it = std::copy(rows.begin(), rows.end(), it)), 42)...};
  return joined_row;
}
} // namespace detail
template <typename... RowViews> auto join(const RowViews&... rows) {
  return detail::join_impl(
      conjunction<
          bool_constant<view_extents_t<RowViews>::rank_dynamic() == 0>...>(),
      rows...);
}

namespace detail {
template <typename Abi, typename V, index E, typename... Vs>
auto load_impl(std::true_type, const Abi&, V, const row_view<E, Vs...>& view,
               index shift) {
  using T = typename V::scalar_type;
  using Vt = typename V::variables_tuple;
  static constexpr std::size_t size = std::tuple_size<Vt>::value;
  std::array<simd<T, Abi>, size> packs;
  int i = 0;
  for_each_tuple_element(
      [&](auto var) {
        packs[i++].copy_from(view[var].data() + shift, element_alignment);
      },
      Vt{});
  return packs;
}

template <typename Abi, typename V, index E, typename... Vs>
auto load_impl(std::false_type, const Abi&, V x, const row_view<E, Vs...>& view,
               index shift) {
  using T = typename variable_traits<V>::value_type;
  simd<T, Abi> pack;
  pack.copy_from(view[x].data() + shift, element_alignment);
  return pack;
}
} // namespace detail

/// \brief Loads all quantities which are refered by a specified view into simd
/// registers.
///
/// Use this to transform a view and position index into a simd vector when
/// doing vecotrised computations.
template <typename Abi, index E, typename... Vs>
add_simd_t<quantities<std::decay_t<Vs>...>, Abi>
load(const row_view<E, Vs...>& view, index shift = 0) noexcept {
  return add_simd_t<quantities<std::decay_t<Vs>...>, Abi>{
      detail::load_impl(is_vector_variable<std::decay_t<Vs>>{}, Abi{},
                        std::decay_t<Vs>{}, view, shift)...};
}

template <index E, typename... Vs>
void store(const quantities<Vs...>& q, const row_view<E, Vs...>& view,
           index shift) noexcept {
  view(shift) = q;
}

template <typename Abi, index E, typename... Vs>
void store(const add_simd_t<quantities<Vs...>, Abi>& q,
           const row_view<E, Vs...>& view, index shift) noexcept {
  auto store = [=](const simd<double, Abi>& pack, double* pointer) {
    pack.copy_to(pointer + shift, element_alignment);
  };
  fub::for_each_tuple_element(
      [&](auto var) { store(q[var], view[var].data()); },
      flatten_variables(Vs{}...));
}

template <typename V, typename Abi> struct simd_proxy;
template <index N, typename Abi, typename... Vars>
struct simd_proxy<row_view<N, Vars...>, Abi> {
  using Variables = std::tuple<Vars...>;
  using simd_type = add_simd_t<quantities<Vars...>, Abi>;

  row_view<N, Vars...> view;
  index shift;

  simd_proxy(row_view<N, Vars...> v, index s) noexcept : view{v}, shift{s} {}
  simd_proxy(const simd_proxy&) /* = delete */;
  simd_proxy& operator=(const simd_proxy&) /* = delete; */;

  operator simd_type() const noexcept { return load<Abi>(view, shift); }

  void operator=(const simd_type& quantities) noexcept {
    store<Abi>(quantities, view, shift);
  }
};

template <typename Abi, index E, typename... Vs>
std::enable_if_t<conjunction<std::is_const<Vs>...>::value,
                 quantities<add_simd_t<std::remove_const_t<Vs>, Abi>...>>
load_or_proxy(const row_view<E, Vs...>& view, index shift = 0) noexcept {
  return load<Abi>(view, shift);
}

template <typename Abi, index E, typename... Vs>
std::enable_if_t<conjunction<negation<std::is_const<Vs>>...>::value,
                 simd_proxy<row_view<E, Vs...>, Abi>>
load_or_proxy(const row_view<E, Vs...>& view, index shift = 0) noexcept {
  return simd_proxy<row_view<E, Vs...>, Abi>{view, shift};
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

template <index Rest, typename F, typename Row, typename... Rows>
F for_each_simd_impl(index_constant<Rest>, F f, const Row& row,
                     const Rows&... rows) {
  return for_each_simd_impl(simd_abi::scalar(), std::move(f), rtake<Rest>(row),
                            rtake<Rest>(rows)...);
}

template <typename F, typename Row, typename... Rows>
F for_each_simd_impl(index_constant<dynamic_extent>, index rest, F f,
                     const Row& row, const Rows&... rows) {
  if (rest > 0) {
    return for_each_simd_impl(simd_abi::scalar(), std::move(f),
                              rtake(row, rest), rtake(rows, rest)...);
  }
  return f;
}

template <typename Abi, typename F, typename FirstRow, typename... Rows>
std::enable_if_t<is_simd_abi_v<Abi>, F>
for_each_simd_impl(std::false_type, const Abi& abi, F f, const FirstRow& first,
                   const Rows&... rows) {
  const index width = simd<double, Abi>::size();
  const index size = first.extents().size();
  const index simd_size = size / width;
  for (index simd_index = 0; simd_index < simd_size; ++simd_index) {
    const index shift = simd_index * width;
    fub::invoke(std::ref(f), abi, load_or_proxy<Abi>(first, shift),
                load_or_proxy<Abi>(rows, shift)...);
  }
  const index rest = size % width;
  return for_each_simd_impl(index_c<dynamic_extent>, rest, std::move(f), first,
                            rows...);
}

template <typename Abi, typename F, typename FirstRow, typename... Rows>
std::enable_if_t<is_simd_abi_v<Abi>, F>
for_each_simd_impl(std::true_type, const Abi& abi, F f, const FirstRow& first,
                   const Rows&... rows) {
  using Extents = view_extents_t<FirstRow>;
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
  static_assert(Extents::rank() == 1,
                "Only one-dimensional views are supported.");
  return for_each_simd_impl(bool_c<Extents::rank_dynamic() == 0>, abi,
                            std::move(f), first, rows...);
}

} // namespace detail

template <typename Abi, typename F, typename... Rows>
std::enable_if_t<conjunction<is_simd_abi<Abi>, is_view<Rows>...>::value, F>
for_each_simd(const Abi& abi, F f, const Rows&... rows) {
  return detail::for_each_simd_impl(abi, std::move(f), rows...);
}

template <typename F, typename... Rows>
std::enable_if_t<conjunction<is_view<Rows>...>::value, F>
for_each_simd(F f, const Rows&... rows) {
  return for_each_simd(simd_abi::native<double>(), std::move(f), rows...);
}

////////////////////////////////////////////////////////////////////////////////
//                                                          [function.permutate]

template <int A, int B, index... Es>
constexpr std::array<index, sizeof...(Es)>
make_permutated_extent_array(const extents<Es...>&) {
  index e[sizeof...(Es)]{Es...};
  index eA = e[A];
  e[A] = e[B];
  e[B] = eA;
  return as_std_array(e);
}

template <int A, int B, typename E, std::size_t... Is>
constexpr auto permutate_extents_impl(const E&, std::index_sequence<Is...>) {
  constexpr std::array<index, sizeof...(Is)> es =
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
template <typename T, index... Es> struct tensor_map<mdspan<T, Es...>> {
  using type = Eigen::TensorMap<Eigen::Tensor<T, sizeof...(Es)>>;
};
template <typename T, index... Es, std::size_t... Is>
tensor_map_t<mdspan<T, Es...>> make_tensor_impl(std::index_sequence<Is...>,
                                                const mdspan<T, Es...>& view) {
  return tensor_map_t<mdspan<T, Es...>>(view.data(), view.extent(Is)...);
}

template <typename T, index... Es>
tensor_map_t<mdspan<T, Es...>> make_tensor(const mdspan<T, Es...>& view) {
  return make_tensor_impl(std::make_index_sequence<sizeof...(Es)>(), view);
}

template <int A, int B, int Rank>
constexpr std::array<int, Rank> make_permutation_index() noexcept {
  auto index = fub::apply(
      [](auto... is) { return std::array<int, sizeof...(is)>{{is()...}}; },
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
        std::array<index, view_rank_v<ViewA>> origin{};
        origin[Dim] = source.extents().get(Dim) - Width;
        auto extents =
            static_cast<std::array<index, view_rank_v<ViewA>>>(dest.extents());
        dt = st.slice(origin, extents);
      },
      view_variables_t<ViewA>());
}

template <int Dim, int Width, typename Extents, typename... Vars>
auto slice_left(const patch_view<Extents, Vars...>& source) {
  auto changed_extents =
      replace_extent(source.extents(), int_c<Dim>, int_c<Width>);
  auto dest =
      make_patch(std::tuple<std::remove_const_t<Vars>...>(), changed_extents);
  slice_left<Dim, Width>(make_view(dest), source);
  return dest;
}

template <int Dim, int Width, typename ViewA, typename ViewB>
void slice_right(const ViewA& dest, const ViewB& source) {
  for_each_tuple_element(
      [&](auto variable) {
        auto st = make_tensor(source[variable]);
        auto dt = make_tensor(dest[variable]);
        std::array<index, view_rank_v<ViewA>> origin{};
        auto extents =
            static_cast<std::array<index, view_rank_v<ViewA>>>(dest.extents());
        dt = st.slice(origin, extents);
      },
      view_variables_t<ViewA>());
}

template <int Dim, int Width, typename Extents, typename... Vars>
auto slice_right(const patch_view<Extents, Vars...>& source) {
  auto changed_extents =
      replace_extent(source.extents(), int_c<Dim>, int_c<Width>);
  auto dest =
      make_patch(std::tuple<std::remove_const_t<Vars>...>(), changed_extents);
  slice_right<Dim, Width>(make_view(dest), source);
  return dest;
}

} // namespace fub

#endif // !BLOCKVIEW_HPP
