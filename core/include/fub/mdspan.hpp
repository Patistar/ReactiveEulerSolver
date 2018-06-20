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

#ifndef FUB_CORE_MDSPAN_HPP
#define FUB_CORE_MDSPAN_HPP

#include "fub/algorithm.hpp"
#include "fub/array.hpp"
#include "fub/extents.hpp"
#include "fub/layout_left.hpp"
#include "fub/optional.hpp"
#include "fub/tuple.hpp"
#include "fub/type_traits.hpp"

#include <array>
#include <cassert>
#include <tuple>
#include <vector>

namespace fub {

template <typename ElementType, typename Extents = extents<dyn>, typename...>
class mdspan;

template <typename T> struct is_mdspan : std::false_type {};
template <typename... Ts> struct is_mdspan<mdspan<Ts...>> : std::true_type {};
template <typename T> static constexpr bool is_mdspan_v = is_mdspan<T>::value;

template <typename T> struct mdspan_extents;
template <typename T, typename E, typename... Ps>
struct mdspan_extents<mdspan<T, E, Ps...>> {
  using type = E;
};
template <typename T> using mdspan_extents_t = typename mdspan_extents<T>::type;

namespace detail {
template <typename ElementT, typename MappingT>
class mdspan_storage : private MappingT {
public:
  using element_type = ElementT;
  using layout_mapping = MappingT;
  using pointer = element_type*;
  using reference = element_type&;

private:
  pointer m_pointer{nullptr};

public:
  constexpr mdspan_storage() = default;

  template <typename... Args>
  constexpr mdspan_storage(pointer ptr, Args&&... args)
      : MappingT{std::forward<Args>(args)...}, m_pointer{ptr} {}

  constexpr const layout_mapping& get_layout_mapping() const noexcept {
    return *this;
  }
  pointer get_pointer() const noexcept { return m_pointer; }
};
} // namespace detail

template <typename ElementT, typename Extents, typename...> class mdspan {
public:
  using element_type = ElementT;
  using value_type = std::remove_const_t<element_type>;
  using extents_type = Extents;
  using pointer = element_type*;
  using reference = element_type&;
  using iterator = pointer;
  using index_type = index;
  using size_type = index;
  using layout = layout_left;
  using layout_mapping = typename layout::template mapping<extents_type>;

private:
  detail::mdspan_storage<element_type, layout_mapping> m_storage{};

public:
  static constexpr index rank = extents_type::rank;
  static constexpr index rank_dynamic = extents_type::rank_dynamic;

  // CONSTRUCTORS

  constexpr mdspan() = default;
  constexpr mdspan(pointer ptr, const extents_type& extents)
      : m_storage(ptr, extents) {}

  // ACCESSOR

  constexpr pointer data() const noexcept { return m_storage.get_pointer(); }

  constexpr index_type size() const noexcept { return extents().size(); }

  constexpr const layout_mapping& get_layout_mapping() const noexcept {
    return m_storage.get_layout_mapping();
  }
  constexpr const extents_type& extents() const noexcept {
    return get_layout_mapping().extents();
  }

  // ELEMENT ACCESS

  template <typename... IndexTypes,
            typename = std::enable_if_t<(sizeof...(IndexTypes) == rank)>,
            typename = std::enable_if_t<
                conjunction<std::is_integral<IndexTypes>...>::value>>
  constexpr reference operator()(IndexTypes... idx) const noexcept {
    return data()[get_layout_mapping()(static_cast<index>(idx)...)];
  }

  constexpr reference operator()(const std::array<index, rank>& indices) const
      noexcept {
    return data()[fub::apply(get_layout_mapping(), indices)];
  }

  constexpr reference operator[](index i) const noexcept { return data()[i]; }

  constexpr operator bool() const noexcept {
    return data() != nullptr;
  }
};

template <typename Mapping, typename Function>
constexpr Function for_each_index(const Mapping& mapping, Function f) {
  constexpr int rank = Mapping::rank;
  using Indices = array<index, rank>;
  optional<Indices> indices = Indices{};
  while (indices) {
    fub::invoke(std::ref(f), *indices);
    indices = mapping.next(*indices);
  }
  return f;
}

/// @brief Invokes the specified function `f` for all indices in the range of
/// the specified extents.
template <typename Function, index... Values>
constexpr Function for_each_index(const extents<Values...>& ext, Function f) {
  using Mapping = layout_left::mapping<extents<Values...>>;
  Mapping mapping{ext};
  return for_each_index(mapping, std::move(f));
}

} // namespace fub

#endif // !SPAN_HPP
