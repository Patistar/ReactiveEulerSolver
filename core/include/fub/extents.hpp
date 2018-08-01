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

/// @file This file introduces `extents<E0, ..., En>`, a compact
/// multi-dimensional size type. Each integral extent `Ei` stands an upper bound
/// in dimension `i` and can either be a compile-time constant signed integral
/// value or `dynamic_extent`. Each compile-time sized extent does not take any
/// extra byte.

#ifndef FUB_CORE_EXTENTS_HPP
#define FUB_CORE_EXTENTS_HPP

#include "fub/dynamic_extent.hpp"
#include "fub/type_traits.hpp"
#include <array>

namespace fub {
inline namespace v1 {
namespace detail {
/// This is the storage type for the extents class and only takes storage for
/// dynamic extents.
///
/// This class is also responsible for accessing all extents.
template <std::size_t Rank, std::size_t DynamicRank,
          std::ptrdiff_t... StaticExtents>
class extents_storage;

/// This is the class template specialisation if all extents are statically
/// known.
///
/// In this case no extra storage is going to be used.
template <std::size_t Rank, std::ptrdiff_t... StaticExtents>
class extents_storage<Rank, 0, StaticExtents...> {
public:
  static_assert(Rank == sizeof...(StaticExtents),
                "Rank does not match sizeof...(StaticExtents)!");
  static_assert(
      conjunction<bool_constant<(StaticExtents != dynamic_extent)>...>::value,
      "Some static extent is equal to dynamic_extent!");

  // [mdspan.extents.obs], Observers of the domain multi-index space

  /// Returns sizeof...(StaticExtents)
  static constexpr std::size_t rank() noexcept { return Rank; }

  /// Returns 0.
  static constexpr std::size_t rank_dynamic() noexcept { return 0; }

  /// Returns the n-th Static Extent.
  static constexpr std::ptrdiff_t static_extent(std::size_t n) noexcept {
    const std::ptrdiff_t static_extents[Rank]{StaticExtents...};
    return static_extents[n];
  }

  /// Returns the N-th Extent.
  ///
  /// In this case where all extents are known at compile time, this returns
  /// `static_extent(n)`.
  constexpr std::ptrdiff_t extent(size_t n) const noexcept {
    return static_extent(n);
  }
};

template <std::size_t Rank, std::ptrdiff_t... StaticExtents>
class extents_storage<Rank, Rank, StaticExtents...> {
public:
  static_assert(Rank == sizeof...(StaticExtents),
                "Rank does not match sizeof...(StaticExtents)!");
  static_assert(
      conjunction<bool_constant<StaticExtents == dynamic_extent>...>::value,
      "Not all static extents are equal to dynamic_extent!");

  // [mdspan.extents.constructors]

  extents_storage() = default;

  template <typename... IndexType,
            typename = std::enable_if_t<conjunction<
                std::is_convertible<IndexType, std::ptrdiff_t>...>::value>,
            typename = std::enable_if_t<sizeof...(IndexType) == Rank>>
  constexpr explicit extents_storage(IndexType... extent) noexcept
      : m_dynamic_extents{static_cast<std::ptrdiff_t>(extent)...} {}

  constexpr explicit extents_storage(
      const std::array<std::ptrdiff_t, Rank>& extents) noexcept
      : m_dynamic_extents{extents} {}

  // [mdspan.extents.obs], Observers of the domain multi-index space

  /// Returns sizeof...(StaticExtents)
  static constexpr std::size_t rank() noexcept { return Rank; }

  /// Returns sizeof...(StaticExtents)
  static constexpr std::size_t rank_dynamic() noexcept { return Rank; }

  /// Returns dynamic_extent.
  static constexpr std::ptrdiff_t static_extent(std::size_t n) noexcept {
    return dynamic_extent;
  }

  /// Returns the N-th Extent.
  constexpr std::ptrdiff_t extent(size_t n) const noexcept {
    return m_dynamic_extents[n];
  }

private:
  std::array<std::ptrdiff_t, Rank> m_dynamic_extents{};
};

template <std::size_t Rank, std::size_t RankDynamic,
          std::ptrdiff_t... StaticExtents>
class extents_storage {
public:
  // [mdspan.extents.constructors]

  extents_storage() = default;

  template <typename... IndexType,
            typename = std::enable_if_t<conjunction<
                std::is_convertible<IndexType, std::ptrdiff_t>...>::value>,
            typename = std::enable_if_t<sizeof...(IndexType) == RankDynamic>>
  constexpr explicit extents_storage(IndexType... extent) noexcept
      : m_dynamic_extents{static_cast<std::ptrdiff_t>(extent)...} {}

  constexpr explicit extents_storage(
      const std::array<std::ptrdiff_t, RankDynamic>& extents) noexcept
      : m_dynamic_extents{extents} {}

  // [mdspan.extents.obs], Observers of the domain multi-index space

  /// Returns sizeof...(StaticExtents)
  static constexpr std::size_t rank() noexcept { return Rank; }

  /// Returns sizeof...(StaticExtents)
  static constexpr std::size_t rank_dynamic() noexcept { return RankDynamic; }

  /// Returns dynamic_extent.
  static constexpr std::ptrdiff_t static_extent(std::size_t n) noexcept {
    constexpr std::ptrdiff_t static_extents[Rank]{StaticExtents...};
    return static_extents[n];
  }

public:
  /// Returns the N-th Extent.
  constexpr std::ptrdiff_t extent(size_t n) const noexcept {
    constexpr std::ptrdiff_t static_extents[Rank]{StaticExtents...};
    if (static_extents[n] != dynamic_extent) {
      return static_extents[n];
    }
    return m_dynamic_extents[find_dynamic_extent_index(n)];
  }

private:
  std::array<std::ptrdiff_t, RankDynamic> m_dynamic_extents{};

  static constexpr std::size_t
  find_dynamic_extent_index(std::size_t n) noexcept {
    constexpr std::ptrdiff_t static_extents[Rank]{StaticExtents...};
    std::size_t dynamic_extent_index = 0;
    for (std::size_t dim = 0; dim < n; ++dim) {
      if (static_extents[dim] == dynamic_extent) {
        ++dynamic_extent_index;
      }
    }
    return dynamic_extent_index;
  }
};

template <typename... IndexType>
constexpr std::size_t count_dynamic_extents(IndexType... extent) noexcept {
  constexpr std::size_t Rank = sizeof...(IndexType);
  std::ptrdiff_t extents[Rank]{static_cast<std::ptrdiff_t>(extent)...};
  std::size_t counter = 0;
  for (std::size_t dim = 0; dim < Rank; ++dim) {
    if (extents[dim] == dynamic_extent) {
      ++counter;
    }
  }
  return counter;
}
} // namespace detail

/// An extents object defines a multidimensional index space which is the
/// Cartesian product of integers extents `[0..N0) * [0..N1) * ...`
template <std::ptrdiff_t... StaticExtents>
class extents
    : private detail::extents_storage<
          sizeof...(StaticExtents),
          detail::count_dynamic_extents(StaticExtents...), StaticExtents...> {
private:
  using base_type =
      detail::extents_storage<sizeof...(StaticExtents),
                              detail::count_dynamic_extents(StaticExtents...),
                              StaticExtents...>;

public:
  // Type Definitions

  using index_type = std::ptrdiff_t;

  // [mdspan.extents.constructors]

  using base_type::base_type;

  // [mdspan.extents.static_observers]

  /// Returns sizeof...(StaticExtents)
  static constexpr std::size_t rank() noexcept { return base_type::rank(); }

  /// Returns sizeof...(StaticExtents)
  static constexpr std::size_t rank_dynamic() noexcept {
    return base_type::rank_dynamic();
  }

  /// Returns the `n`-th `StaticExtent`.
  static constexpr std::ptrdiff_t static_extent(std::size_t n) noexcept {
    return base_type::static_extent(n);
  }

  // [mdspan.extents.observers]

  /// Returns the `n`-th run-time extent.
  constexpr std::ptrdiff_t extent(size_t n) const noexcept {
    return base_type::extent(n);
  }
};

// [mdspan.extents.traits.is_extent]

/// \ingroup type-traits
/// This is true `std::true_type` iff `E` is `extents<Es...>`
/// for some `std::ptrdiff_t... Es`.
/// @{
template <typename E> struct is_extents : std::false_type {};
template <std::ptrdiff_t... StaticExtents>
struct is_extents<extents<StaticExtents...>> : std::true_type {};
#ifdef FUB_WITH_INLINE_VARIABLES
template <typename E> inline constexpr bool is_extents_v = is_extents<E>::value;
#else
template <typename E> static constexpr bool is_extents_v = is_extents<E>::value;
#endif
/// @}

/// Returns: `true` if `lhs.rank() == rhs.rank()` and `lhs.extents(r) ==
/// rhs.extents(r)` for all `r` in the range `[0, lhs.rank())`, or false
/// otherwise.
template <std::ptrdiff_t... StaticExtentsL, std::ptrdiff_t... StaticExtentsR>
constexpr bool operator==(const extents<StaticExtentsL...>& lhs,
                          const extents<StaticExtentsR...>& rhs) noexcept {
  if (lhs.rank() != rhs.rank()) {
    return false;
  }
  for (std::size_t r = 0; r < lhs.rank(); ++r) {
    if (lhs.extent(r) != rhs.extent(r)) {
      return false;
    }
  }
  return true;
}

template <std::ptrdiff_t... StaticExtentsL, std::ptrdiff_t... StaticExtentsR>
constexpr bool operator!=(const extents<StaticExtentsL...>& left,
                          const extents<StaticExtentsR...>& right) noexcept {
  return !(left == right);
}

template <std::ptrdiff_t... StaticExtents>
constexpr std::ptrdiff_t size(const extents<StaticExtents...>& e) noexcept {
  std::ptrdiff_t count = 1;
  for (std::size_t dim = 0; dim < e.rank(); ++dim) {
    count *= e.extent(dim);
  }
  return count;
}

namespace detail {
template <std::ptrdiff_t... StaticExtents, std::size_t... Is>
constexpr std::array<std::ptrdiff_t, sizeof...(StaticExtents)>
as_array_impl(const extents<StaticExtents...>& e,
              std::index_sequence<Is...>) noexcept {
  return {{e.extent(Is)...}};
}

template <typename T, std::size_t N, std::size_t... Is>
constexpr std::array<T, N>
as_std_array_impl(const T (&array)[N], std::index_sequence<Is...>) noexcept {
  return {{array[Is]...}};
}
} // namespace detail

template <typename Extents>
using index_array_t = std::array<typename Extents::index_type, Extents::rank()>;

template <typename Extents>
constexpr bool is_in_range(const Extents& e,
                           const index_array_t<Extents>& array) noexcept {
  for (std::size_t r = 0; r < Extents::rank(); ++r) {
    if (array[r] < 0 || e.extent(r) <= array[r]) {
      return false;
    }
  }
  return true;
}

template <std::ptrdiff_t... StaticExtents>
constexpr std::array<std::ptrdiff_t, sizeof...(StaticExtents)>
as_array(const extents<StaticExtents...>& e) noexcept {
  return detail::as_array_impl(
      e, std::make_index_sequence<sizeof...(StaticExtents)>());
}

template <typename T, std::size_t N>
constexpr std::array<T, N> as_std_array(const T (&array)[N]) noexcept {
  return detail::as_std_array_impl(array, std::make_index_sequence<N>());
}

template <typename E,
          typename std::enable_if_t<is_extents<E>::value, int*> = nullptr>
constexpr std::array<std::ptrdiff_t, E::rank_dynamic()> filter_dynamic_extents(
    const std::array<std::ptrdiff_t, E::rank()>& array) noexcept {
  constexpr std::size_t RankDynamic = E::rank_dynamic();
  std::ptrdiff_t dynamic_extents[RankDynamic]{};
  std::size_t counter = 0;
  for (std::size_t r = 0; r < E::rank(); ++r) {
    if (E::static_extent(r) == dynamic_extent) {
      dynamic_extents[counter++] = array[r];
    }
  }
  assert(counter == RankDynamic);
  return as_std_array(dynamic_extents);
}

template <std::ptrdiff_t... StaticExtents>
constexpr std::array<std::ptrdiff_t, extents<StaticExtents...>::rank_dynamic()>
get_dynamic_extents(const extents<StaticExtents...>& e) noexcept {
  return filter_dynamic_extents<extents<StaticExtents...>>(as_array(e));
}

////////////////////////////////////////////////////////////////////////////////
// grow
// {{{
namespace detail {
template <std::size_t I, std::size_t J, index Extent, index Width>
struct grow_extent {
  static constexpr index value = Extent;
};

template <std::size_t I, index Extent, index Width>
struct grow_extent<I, I, Extent, Width> {
  static_assert(Extent > 0 && Width >= 0,
                "This is only allowed to be instantiated if Extent and Width "
                "are positive.");
  static constexpr index value = Extent + Width;
};

template <std::ptrdiff_t... StaticExtents> struct static_grow_impl {
  using Result = extents<StaticExtents...>;
  template <typename E,
            typename std::enable_if_t<(E::rank_dynamic() > 0), int*> = nullptr>
  constexpr Result operator()(const E& e) noexcept {
    return Result(get_dynamic_extents(e));
  }

  template <typename E,
            typename std::enable_if_t<(E::rank_dynamic() == 0), int*> = nullptr>
  constexpr Result operator()(const E& e) noexcept {
    return Result();
  }
};

template <int Dim, index... Es, std::size_t... Is>
constexpr auto grow_impl(const extents<Es...>& e, int_constant<Dim>,
                         std::false_type, std::index_sequence<Is...>) noexcept {
  return static_grow_impl<grow_extent<Is, Dim, Es, 1>::value...>{}(e);
}

/// Grow for the case that Dim is a static extent
template <int Dim, index... Es>
constexpr auto grow_impl(const extents<Es...>& e, int_constant<Dim> dim,
                         std::false_type) noexcept {
  constexpr std::size_t Rank = extents<Es...>::rank();
  return grow_impl(e, dim, std::false_type{}, std::make_index_sequence<Rank>());
}

template <int Dim, index... Es, std::size_t... Is>
constexpr extents<Es...> grow_impl(const extents<Es...>& e, int_constant<Dim>,
                                   std::true_type,
                                   std::index_sequence<Is...>) noexcept {
  using E = extents<Es...>;
  index native[E::rank()]{e.extent(Is)...};
  native[Dim] += 1;
  const std::array<index, E::rank_dynamic()> dynamic_extents =
      filter_dynamic_extents<E>(as_std_array(native));
  return E(dynamic_extents);
}

template <int Dim, index... Es>
constexpr extents<Es...> grow_impl(const extents<Es...>& e, int_constant<Dim>,
                                   std::true_type) noexcept {
  constexpr std::size_t rank = extents<Es...>::rank();
  return grow_impl(e, int_c<Dim>, std::true_type{},
                   std::make_index_sequence<rank>{});
}
} // namespace detail

template <int Dim, index... Es>
constexpr auto grow(const extents<Es...>& e,
                    int_constant<Dim> dim = int_constant<Dim>()) noexcept {
  return detail::grow_impl(
      e, dim, bool_c<extents<Es...>::static_extent(Dim) == dynamic_extent>);
}
// }}}

////////////////////////////////////////////////////////////////////////////////
// replace_extent
// {{{
namespace detail {
template <std::size_t ReplaceIndex, index ReplaceWith, std::size_t Visitor,
          index Value>
struct replace_index : index_constant<Value> {};

template <std::size_t ReplaceIndex, index ReplaceWith, index Value>
struct replace_index<ReplaceIndex, ReplaceWith, ReplaceIndex, Value>
    : index_constant<ReplaceWith> {};

template <std::ptrdiff_t... StaticExtents> struct replace_extent_impl {
  using Result = extents<StaticExtents...>;

  template <typename E>
  constexpr Result operator()(const E&, size_constant<0>) {
    return Result{};
  }

  template <typename E, std::size_t RankDynamic>
  constexpr Result operator()(const E& extents, size_constant<RankDynamic>) {
    static_assert(Result::rank() == E::rank(), "");
    static_assert(Result::rank_dynamic() == RankDynamic, "");
    static_assert(Result::rank_dynamic() <= E::rank_dynamic(), "");
    return Result{filter_dynamic_extents<Result>(as_array(extents))};
  }
};

template <int Dim, int Value, index... Es, std::size_t... Is>
constexpr auto replace_extent(const extents<Es...>& e, int_constant<Dim>,
                              int_constant<Value>,
                              std::index_sequence<Is...>) noexcept {
  using Result = extents<replace_index<Dim, Value, Is, Es>::value...>;
  return replace_extent_impl<replace_index<Dim, Value, Is, Es>::value...>{}(
      e, size_c<Result::rank_dynamic()>);
}
} // namespace detail

template <int Dim, int Value, index... Es>
constexpr auto
replace_extent(const extents<Es...>& e,
               int_constant<Dim> dim = int_constant<Dim>(),
               int_constant<Value> value = int_constant<Value>()) noexcept {
  return detail::replace_extent(e, dim, value,
                                std::make_index_sequence<sizeof...(Es)>());
}

template <typename E, int Dim, int V>
using replace_extent_t =
    decltype(replace_extent(std::declval<E>(), int_c<Dim>, int_c<V>));
// }}}

namespace detail {
template <int Rank> struct dynamic_extents {
  template <typename> struct impl;
  template <int... Is> struct impl<std::integer_sequence<int, Is...>> {
    using type = extents<(Is, dynamic_extent)...>;
  };
  using type = typename impl<std::make_integer_sequence<int, Rank>>::type;
};
} // namespace detail
template <int Rank>
using dynamic_extents_t = typename detail::dynamic_extents<Rank>::type;

} // namespace v1
} // namespace fub

#endif
