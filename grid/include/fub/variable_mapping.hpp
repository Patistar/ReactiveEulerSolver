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

#ifndef FUB_GRID_VARIABLE_MAPPING_HPP
#define FUB_GRID_VARIABLE_MAPPING_HPP

#include "fub/variable_list.hpp"
#include "fub/variable_ref.hpp"

namespace fub {
inline namespace v1 {
template <typename M> struct mdspan_add_const;

template <typename T, typename E, typename L, typename A>
struct mdspan_add_const<basic_mdspan<T, E, L, A>> {
  using type = basic_mdspan<std::add_const_t<T>, E, L,
                            typename A::template rebind<std::add_const_t<T>>>;
};

template <typename M>
using mdspan_add_const_t = typename mdspan_add_const<M>::type;

template <typename M> struct mdspan_remove_const;
template <typename T, typename E, typename L, typename A>
struct mdspan_remove_const<basic_mdspan<T, E, L, A>> {
  using type =
      basic_mdspan<std::remove_const_t<T>, E, L,
                   typename A::template rebind<std::remove_const_t<T>>>;
};
template <typename M>
using mdspan_remove_const_t = typename mdspan_remove_const<M>::type;

/// This is a utility class which helps to map multidimensional variable data
/// in some congtiguous range of memory.
template <typename VariableList, typename MdSpan> struct variable_mapping {
  using mdspan_type = MdSpan;
  using const_mdspan_type = mdspan_add_const_t<MdSpan>;
  using value_type = typename mdspan_type::value_type;
  using const_value_type = std::add_const_t<value_type>;
  using mapping_type = typename mdspan_type::mapping;
  using extents_type = typename mdspan_type::extents_type;
  using accessor_type = rebind_t<typename mdspan_type::accessor, value_type>;
  using const_accessor_type = rebind_t<accessor_type, const_value_type>;
  using layout = typename mdspan_type::layout;

  static constexpr std::ptrdiff_t static_size() noexcept {
    constexpr std::ptrdiff_t vsize = VariableList::static_size();
    constexpr std::ptrdiff_t ssize =
        static_required_span_size(layout(), extents_type());
    if (vsize == dynamic_extent || ssize == dynamic_extent) {
      return dynamic_extent;
    }
    return vsize * ssize;
  }

  using span_type = basic_span<value_type, static_size(), accessor_type>;
  using const_span_type = typename const_mdspan_type::span_type;

  template <typename Tag>
  static constexpr mdspan_type make_mdspan_(std::true_type, span_type s,
                                            extents_type e, VariableList vars,
                                            Tag tag) {
    optional<std::ptrdiff_t> index = vars.index(tag);
    assert(index);
    const mapping_type mapping(e);
    const std::ptrdiff_t size = mapping.required_span_size();
    const std::ptrdiff_t lower = (*index) * size;
    return mdspan_type(subspan(s, lower, size), mapping);
  }

  template <typename Tag>
  static constexpr const_mdspan_type
  make_mdspan_(std::false_type, const_span_type s, extents_type e,
               VariableList vars, Tag tag) {
    optional<std::ptrdiff_t> index = vars.index(tag);
    assert(index);
    const mapping_type mapping(e);
    const std::ptrdiff_t size = mapping.required_span_size();
    const std::ptrdiff_t lower = (*index) * size;
    return const_mdspan_type(subspan(s, lower, size), mapping);
  }

  template <typename SpanLike, typename Tag,
            typename = std::enable_if_t<
                std::is_convertible<SpanLike, const_span_type>{}>>
  static constexpr auto make_mdspan(SpanLike s, extents_type e,
                                    VariableList vars, Tag tag) {
    return make_mdspan_(std::is_convertible<SpanLike, span_type>{}, s, e, vars,
                        tag);
  }
};

} // namespace v1
} // namespace fub

#endif