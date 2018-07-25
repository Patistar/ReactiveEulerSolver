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

#include "fub/mdspan.hpp"

namespace fub {
///////////////////////////////////////////////////////////////////////////////
// Forward Declerations

template <typename Vs, typename E, typename L, typename A>
class basic_patch_view;

///////////////////////////////////////////////////////////////////////////////
//                                                           [patch_view.class]

template <typename Span, typename Mapping>
class basic_patch_view_storage : Mapping {
public:
  basic_patch_view_storage() = default;
  basic_patch_view_storage(Span span, Mapping mapping) noexcept
      : Mapping(mapping), m_span(span) {}

  constexpr Span span() const noexcept { return m_span; }

  constexpr Mapping get_mapping() const noexcept { return *this; }

private:
  Span m_span;
};

template <typename Mapping, typename VariableList> struct compute_span_extents {
  static constexpr std::ptrdiff_t
  extent(hana::basic_type<Mapping>, hana::basic_type<VariableList>) noexcept {
    if (Mapping::static_required_span_size() == dynamic_extent ||
        VariableList::static_size() == dynamic_extent) {
      return dynamic_extent;
    }
    return Mapping::static_required_span_size() * VariableList::static_size();
  }

  using type =
      extents<extent(hana::type_c<Mapping>, hana::type_c<VariableList>)>;
};

template <typename Mapping, typename VariableList>
using compute_span_extents_t =
    typename compute_span_extents<Mapping, VariableList>::type;

/// This a view type for `patch<Variables, Extents, Storage>` which provides
/// multi dimensional views for distinct variables with uniform extents, layout
/// and accessor policies.
///
/// *Example:*
/// \code
///     struct MyVariable : scalar_variable {};
///     tag_t<MyVariable> variable;
///
///     // Create a fixed-size patch containing data for variable MyVariable
///     patch<MyVariable, extents<32, 32>> p;
///
///     // Create a view on the patch p.
///     patch_view<MyVariable, extents<32, 32>> view = p;
///
///     // Assign data
///     view[variable](10, 10) = 42.;
///     
///     // Access data
///     REQUIRE(view[variable](10, 10) == 42.);
/// \endcode
template <typename VariableList, typename MdSpan>
class basic_patch_view : VariableList {
public:
  using variable_list = VariableList;
  using extents = Extents;
  using accessor = Accessor;
  using layout_policy = LayoutPolicy;
  using mapping = typename LayoutPolicy::template mapping<Extents>;
  using element_type = typename Accessor::element_type;
  using pointer = typename Accessor::pointer;
  using mdspan = basic_mdspan<element_type, extents, layout_policy, accessor>;
  using span_extents = compute_span_extents_t<mapping, variable_list>;
  using span_type =
      basic_span<element_type, span_extents::static_extent(0), accessor>;

  template <typename Tag>
  static constexpr bool is_valid_tag_v =
      variable_list::template is_valid_tag_v<Tag>;

  basic_patch_view() = default;

  constexpr basic_patch_view(span_type data, variable_list list,
                             extents e = extents()) noexcept
      : VariableList(list), m_storage(data, mapping(e)) {
    assert(span().size() >= required_span_size());
  }

  template <typename Variable, std::ptrdiff_t N,
            typename std::enable_if_t<is_valid_tag_v<Variable>>* = nullptr>
  constexpr mdspan operator[](basic_tag<Variable, N> tag) const noexcept {
    const variable_list variables = get_variable_list();
    optional<std::ptrdiff_t> index = variables.index(tag);
    assert(index);
    mapping map = get_mapping();
    const std::ptrdiff_t size = map.required_span_size();
    return mdspan(subspan(span(), *index * size, size), map);
  }

  constexpr variable_list get_variable_list() const noexcept { return *this; }

  constexpr span_type span() const noexcept { return m_storage.span(); }

  constexpr mapping get_mapping() const noexcept {
    return m_storage.get_mapping();
  }

  constexpr std::ptrdiff_t required_span_size() const noexcept {
    return get_variable_list().size() * get_mapping().required_span_size();
  }

  friend constexpr bool operator==(const basic_patch_view& lhs,
                                   const basic_patch_view& rhs) noexcept {
    return lhs.span() == rhs.span() &&
           lhs.get_variable_list() == rhs.get_variable_list() &&
           lhs.get_mapping() == rhs.get_mapping();
  }

  friend constexpr bool operator!=(const basic_patch_view& lhs,
                                   const basic_patch_view& rhs) noexcept {
    return !(lhs == rhs);
  }

private:
  basic_patch_view_storage<span_type, mapping> m_storage;
};

template <typename VariableList, typename T, typename Extents>
using variable_mdspan =
    basic_patch_view<VariableList, T, Extents, layout_left, accessor_native<T>>;

} // namespace fub

#endif // !BLOCKVIEW_HPP
