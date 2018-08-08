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

#include "fub/accessor_native.hpp"
#include "fub/extents.hpp"
#include "fub/layout_left.hpp"
#include "fub/span.hpp"
#include "fub/tuple.hpp"
#include "fub/type_traits.hpp"

#include <array>
#include <cassert>
#include <tuple>
#include <vector>

namespace fub {
inline namespace v1 {

// Forward Decleration of `mdspan`
template <typename T, typename E, typename L, typename A> class basic_mdspan;

///////////////////////////////////////////////////////////////////////////////
//                                                     [mdspan.trait.is_mdspan]

template <typename T> struct is_mdspan : std::false_type {};
template <typename T, typename E, typename L, typename A>
struct is_mdspan<basic_mdspan<T, E, L, A>> : std::true_type {};
template <typename T> static constexpr bool is_mdspan_v = is_mdspan<T>::value;

///////////////////////////////////////////////////////////////////////////////
//                                                [mdspan.trait.mdspan_extents]

template <typename T> struct mdspan_extents;
template <typename T, typename E, typename L, typename A>
struct mdspan_extents<basic_mdspan<T, E, L, A>> {
  using type = E;
};
template <typename T> using mdspan_extents_t = typename mdspan_extents<T>::type;

///////////////////////////////////////////////////////////////////////////////
//                                                             [mdspan.storage]

namespace detail {
template <typename T, typename Mapping, typename Accessor>
class basic_mdspan_storage : Mapping, Accessor {
public:
  using element_type = T;
  using mapping = Mapping;
  using accessor = Accessor;
  using pointer = typename accessor::pointer;
  using reference = typename accessor::reference;

private:
  pointer m_pointer{};

public:
  constexpr basic_mdspan_storage() = default;

  template <typename... Args>
  constexpr explicit basic_mdspan_storage(pointer ptr, mapping m = mapping(),
                                          accessor a = accessor())
      : mapping(m), accessor(a), m_pointer(ptr) {}

  constexpr const mapping& get_mapping() const noexcept { return *this; }

  constexpr const accessor& get_accessor() const noexcept { return *this; }

  pointer get_pointer() const noexcept { return m_pointer; }

  reference access(std::ptrdiff_t n) const noexcept {
    return get_accessor().access(m_pointer, get_mapping().required_span_size(),
                                 n);
  }
};
} // namespace detail

///////////////////////////////////////////////////////////////////////////////
//                                                               [mdspan.class]

template <typename T, typename Extents, typename Layout = layout_left,
          typename Accessor = accessor_native<T>>
class basic_mdspan {
public:
  using element_type = T;
  using extents_type = Extents;
  using layout = Layout;
  using accessor = Accessor;
  using value_type = std::remove_const_t<element_type>;
  using pointer = typename accessor::pointer;
  using reference = typename accessor::reference;
  using index_type = typename Extents::index_type;
  using size_type = std::size_t;
  using mapping = typename layout::template mapping<Extents>;
  using span_type =
      basic_span<element_type,
                 static_required_span_size(layout(), extents_type()), accessor>;

  // CONSTRUCTORS

  constexpr basic_mdspan() = default;
  basic_mdspan(const basic_mdspan&) = default;

  template <
      typename... Args,
      typename std::enable_if_t<std::is_constructible<Extents, Args...>::value,
                                void*> = nullptr>
  constexpr basic_mdspan(pointer ptr, Args&&... args)
      : m_storage(ptr, mapping(extents_type(std::forward<Args>(args)...))) {}

  constexpr basic_mdspan(pointer ptr, const mapping& m,
                         const accessor& a = accessor())
      : m_storage(ptr, m, a) {}

  constexpr basic_mdspan(span_type span, const mapping& m)
      : m_storage(span.data(), m, span.get_accessor()) {}

  template <
      typename S, typename OtherExtents, typename OtherLayout,
      typename OtherAccess,
      typename = std::enable_if_t<std::is_constructible<
          span_type, typename basic_mdspan<S, OtherExtents, OtherLayout,
                                           OtherAccess>::span_type>::value>>
  constexpr basic_mdspan(
      const basic_mdspan<S, OtherExtents, OtherLayout, OtherAccess>& other)
      : basic_mdspan(other.span(), mapping(other.get_extents())) {}


  // [mdspan.basic.domobs], basic_mdspan observers of the domain multi-index
  // space

  static constexpr index_type rank() noexcept { return Extents::rank(); }
  static constexpr index_type rank_dynamic() noexcept {
    return Extents::rank_dynamic();
  }
  static constexpr index_type static_extent(std::size_t r) noexcept {
    return Extents::static_extent();
  }

  constexpr Extents get_extents() const noexcept {
    return get_mapping().get_extents();
  }

  constexpr index_type extent(std::size_t n) const noexcept {
    return get_extents().extent(n);
  }
  constexpr index_type size() const noexcept {
    using fub::size;
    return size(get_extents());
  }
  constexpr index_type unique_size() const noexcept {
    return get_mapping().unique_size();
  }

  // [mdspan.basic.obs], basic_mdspan observers of the mapping

  static constexpr bool is_always_unique() noexcept {
    return mapping::is_always_unique();
  }
  static constexpr bool is_always_contiguous() noexcept {
    return mapping::is_always_contiguous();
  }
  static constexpr bool is_always_strided() noexcept {
    return mapping::is_always_strided();
  }

  constexpr mapping get_mapping() const noexcept {
    return m_storage.get_mapping();
  }
  constexpr bool is_unique() const noexcept {
    return get_mapping().is_unique();
  }
  constexpr bool is_contiguous() const noexcept {
    return get_mapping().is_contiguous();
  }
  constexpr bool is_strided() const noexcept {
    return get_mapping().is_strided();
  }
  constexpr index_type stride(size_t rank) const {
    return get_mapping().stride(rank);
  }

  // [mdspan.basic.codomain], basic_mdspan observers of the codomain

  constexpr accessor get_accessor() const noexcept {
    return m_storage.get_accessor();
  }

  constexpr basic_span<element_type,
                       static_required_span_size(layout(), extents_type()),
                       accessor>
  span() const noexcept {
    return {m_storage.get_pointer(), size(), get_accessor()};
  }

  // ELEMENT ACCESS

  template <typename... IndexTypes,
            typename = std::enable_if_t<(sizeof...(IndexTypes) == rank())>,
            typename = std::enable_if_t<conjunction<
                std::is_convertible<IndexTypes, index_type>...>::value>>
  constexpr reference operator()(IndexTypes... indices) const noexcept {
    return m_storage.access(get_mapping()(indices...));
  }

  constexpr reference operator()(const std::array<index, rank()>& indices) const
      noexcept {
    return m_storage.access(fub::apply(get_mapping(), indices));
  }

  constexpr reference operator[](index i) const noexcept {
    return m_storage.access(i);
  }

  constexpr operator bool() const noexcept {
    return m_storage.get_pointer() != nullptr;
  }

private:
  detail::basic_mdspan_storage<element_type, mapping, accessor> m_storage{};
};

template <typename T, std::ptrdiff_t... Extents>
using mdspan =
    basic_mdspan<T, extents<Extents...>, layout_left, accessor_native<T>>;

} // namespace v1
} // namespace fub

#endif // !SPAN_HPP
