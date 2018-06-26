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
/// value or `dyn`. Each compile-time sized extent does not take any extra byte.
///
/// Example 1:
///     constexpr extents<2, 4> e;
///     static_assert(sizeof(e) == sizeof(char));
///     static_assert(e.get(0) == 2);
///     static_assert(e.get(1) == 4);
///     static_assert(e.size() == 8);
///
/// Example 2:
///     constexpr extents<2, dyn> e(4);
///     static_assert(sizeof(e) == sizeof(std::ptrdiff_t));
///     static_assert(e.get(0) == 2);
///     static_assert(e.get(1) == 4);
///     static_assert(e.size() == 8);
///
/// Example 3:
///     constexpr extents<dyn, dyn> e(2, 4);
///     static_assert(sizeof(e) == 2 * sizeof(std::ptrdiff_t));
///     static_assert(e.get(0) == 2);
///     static_assert(e.get(1) == 4);
///     static_assert(e.size() == 8);

#ifndef FUB_CORE_EXTENTS_HPP
#define FUB_CORE_EXTENTS_HPP

#include "fub/algorithm.hpp"
#include "fub/tuple.hpp"
#include "fub/type_traits.hpp"

#include <range/v3/algorithm/count.hpp>
#include <range/v3/algorithm/equal.hpp>

#include <array>

namespace fub {
template <typename E> struct extents_rank {
  static constexpr int value = E::rank;
};
template <typename E>
static constexpr int extents_rank_v = extents_rank<E>::value;

template <typename E> struct extents_rank_dynamic {
  static constexpr int value = E::rank_dynamic;
};
template <typename E>
static constexpr int extents_rank_dynamic_v = extents_rank_dynamic<E>::value;

namespace detail {

template <int Rank, int Dynamic, index... Es> class extents_;

/// @brief Extents implementation if all extents are compile-time known.
template <int Rank, index... Es> class extents_<Rank, 0, Es...> {
public:
  static constexpr int rank_dynamic = 0;
  static constexpr int rank = Rank;

  extents_() = default;

  /// @brief This Constructor throws an exception if any specified extent is not
  /// equal to the compile time known extent.
  constexpr explicit extents_(const std::array<index, Rank>& extents) {
    std::array<index, rank> check{{Es...}};
    if (!ranges::equal(check, extents)) {
      throw std::invalid_argument{
          "Specified extents do not match the compile time known extents."};
    }
  }

  constexpr explicit operator std::array<index, rank>() const noexcept {
    return {{Es...}};
  }

  constexpr index size() const noexcept {
    constexpr std::array<index, rank> extents{{Es...}};
    return fub::accumulate(extents, index(1), std::multiplies<>());
  }

  static constexpr index static_get(int dim) noexcept {
    constexpr std::array<index, rank> es{{Es...}};
    return es[dim];
  }

  constexpr index get(int dim) const noexcept { return static_get(dim); }

  template <typename Archive> void serialize(Archive&, unsigned) {}
};

/// @brief Extents implementation if all extents are runtime known.
template <int Rank, index... Es> class extents_<Rank, Rank, Es...> {
  std::array<index, Rank> m_extents;

public:
  static constexpr int rank_dynamic = Rank;
  static constexpr int rank = Rank;

  /// @brief Constructs an empty extents of size zero.
  constexpr extents_() = default;

  constexpr explicit extents_(const std::array<index, Rank>& extents)
      : m_extents{extents} {}

  template <typename... Is, typename = std::enable_if_t<
                                conjunction<std::is_integral<Is>...>::value &&
                                (sizeof...(Is) == Rank)>>
  constexpr explicit extents_(Is... extents)
      : m_extents{{static_cast<index>(extents)...}} {}

  constexpr explicit operator std::array<index, Rank>() const noexcept {
    return m_extents;
  }

  constexpr index size() const noexcept {
    return fub::accumulate(m_extents.begin(), m_extents.end(), index(1),
                           std::multiplies<>());
  }

  static constexpr index static_get(int) noexcept { return dyn; }

  constexpr index get(int dim) const noexcept { return m_extents[dim]; }

  template <typename Archive> void serialize(Archive& ar, unsigned) {
    ar& m_extents;
  }
};

/// @brief Extents implementation for the case where we have mixed
/// compile-time AND runtime sizes.
///
/// In this case we store the runtime extents in an internal std::array and map
/// the requested dimension `get(dim)` to the proper internal index.
///
/// @todo Construction of an extent of this type is not canonic.
template <int Rank, int Dynamic, index... Es> class extents_ {
  std::array<index, Dynamic> m_dynamic_extents;

public:
  static constexpr int rank_dynamic = Dynamic;
  static constexpr int rank = Rank;

  static_assert(rank_dynamic < rank);

  constexpr extents_() = default;

  template <typename... Is, typename = std::enable_if_t<
                                conjunction<std::is_integral<Is>...>::value &&
                                (sizeof...(Is) == rank_dynamic)>>
  constexpr extents_(Is... dynamic_extents)
      : m_dynamic_extents{{static_cast<index>(dynamic_extents)...}} {}

  constexpr extents_(const std::array<index, rank>& a) : m_dynamic_extents{} {
    constexpr std::array<index, rank> es{{Es...}};
    int dynamic_extents_counter{0};
    for (int dim = 0; dim < rank; ++dim) {
      if (es[dim] == dyn) {
        m_dynamic_extents[dynamic_extents_counter] = a[dim];
        ++dynamic_extents_counter;
      } else {
        if (es[dim] != a[dim]) {
          throw std::invalid_argument{
              "Specified extents do not match the compile time known extents."};
        }
      }
    }
  }

  static constexpr index static_get(int dim) noexcept {
    std::array<index, rank> es{{Es...}};
    return es[dim];
  }

  constexpr index get(int dim) const noexcept {
    std::array<index, rank> es{{Es...}};
    if (es[dim] == dyn) {
      int count = fub::count(es.begin(), es.begin() + dim, dyn);
      return m_dynamic_extents[count];
    } else {
      return es[dim];
    }
  }

private:
  template <int... Is>
  constexpr auto as_array(std::integer_sequence<int, Is...>) const noexcept {
    return std::array<index, rank>{{get(Is)...}};
  }

public:
  constexpr explicit operator std::array<index, rank>() const noexcept {
    return as_array(make_int_sequence<rank>());
  }

  constexpr index size() const noexcept {
    auto es = static_cast<std::array<index, rank>>(*this);
    return fub::accumulate(es, index(1), std::multiplies<>());
  }

  template <typename Archive> void serialize(Archive& ar, unsigned) {
    ar& m_dynamic_extents;
  }
};

} // namespace detail

template <index... Es>
class extents
    : private detail::extents_<sizeof...(Es), fub::count({Es...}, dyn), Es...> {
public:
  using base = detail::extents_<sizeof...(Es), fub::count({Es...}, dyn), Es...>;
  static constexpr int rank = base::rank;
  static constexpr int rank_dynamic = base::rank_dynamic;

  using base::base;
  using base::get;
  using base::serialize;
  using base::size;
  using base::static_get;

  using base::operator std::array<index, rank>;
};

template <index... Es, index... Fs>
constexpr std::enable_if_t<(sizeof...(Es) == sizeof...(Fs)), bool>
operator==(const extents<Es...>& e1, const extents<Fs...>& e2) noexcept {
  constexpr int rank = sizeof...(Es);
  auto a1 = static_cast<std::array<index, rank>>(e1);
  auto a2 = static_cast<std::array<index, rank>>(e2);
  return ranges::equal(a1, a2);
}

template <index... Es, index... Fs>
constexpr std::enable_if_t<(sizeof...(Es) == sizeof...(Fs)), bool>
operator!=(const extents<Es...>& e1, const extents<Fs...>& e2) noexcept {
  return !(e1 == e2);
}

template <index... Es>
constexpr auto as_array(const extents<Es...>& e) noexcept {
  return static_cast<std::array<index, sizeof...(Es)>>(e);
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
  static constexpr index value = Extent + Width;
};

template <int Dim, index... Es, int... Is>
constexpr auto grow_impl(const extents<Es...>& e, int_constant<Dim>,
                         std::false_type,
                         std::integer_sequence<int, Is...>) noexcept {
  auto a = static_cast<std::array<index, sizeof...(Es)>>(e);
  a[Dim] += 1;
  return extents<grow_extent<Dim, Is, Es, 1>::value...>(a);
}

/// Grow for the case that Dim is a static extent
template <int Dim, index... Es>
constexpr auto grow_impl(const extents<Es...>& e, int_constant<Dim> dim,
                         std::false_type) noexcept {
  constexpr int rank = extents<Es...>::rank;
  return grow_impl(e, dim, std::false_type{}, make_int_sequence<rank>());
}

template <int Dim, index... Es>
constexpr extents<Es...> grow_impl(const extents<Es...>& e, int_constant<Dim>,
                                   std::true_type) noexcept {
  constexpr int rank = extents<Es...>::rank;
  auto a = static_cast<std::array<index, rank>>(e);
  a[Dim] += 1;
  return extents<Es...>(a);
}
} // namespace detail

template <int Dim, index... Es>
constexpr auto grow(const extents<Es...>& e, int_constant<Dim> dim) noexcept {
  return detail::grow_impl(e, dim,
                           bool_c<extents<Es...>::static_get(Dim) == dyn>);
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

template <int Dim, int Value, index... Es, std::size_t... Is>
constexpr auto replace_extent(const extents<Es...>&, int_constant<Dim>,
                              int_constant<Value>,
                              std::index_sequence<Is...>) noexcept {
  std::array<index, sizeof...(Is)> replaced{
      {replace_index<Dim, Value, Is, Es>::value...}};
  return extents<replace_index<Dim, Value, Is, Es>::value...>(replaced);
}
} // namespace detail

template <int Dim, int Value, index... Es>
constexpr auto replace_extent(const extents<Es...>& e, int_constant<Dim> dim,
                              int_constant<Value> value) noexcept {
  return detail::replace_extent(e, dim, value,
                                std::make_index_sequence<sizeof...(Es)>());
}

template <typename E, int Dim, int V>
using replace_extent_t =
    decltype(replace_extent(std::declval<E>(), int_c<Dim>, int_c<V>));
// }}}

} // namespace fub

#endif
