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

#ifndef FUB_GRID_VARIABLE_VIEW_HPP
#define FUB_GRID_VARIABLE_VIEW_HPP

#include "fub/variable_mapping.hpp"

namespace fub {
inline namespace v1 {
/// This is view on contiguous data but indexed by a variable list and a uniform
/// multi-dimensional mapping.
/// \nosubgrouping
template <typename VariableList, typename MdSpan>
class basic_variable_view : VariableList, MdSpan::mapping {
public:
  /// \name Types

  using mdspan_type = MdSpan;
  using value_type = typename mdspan_type::value_type;
  using accessor_type = typename mdspan_type::accessor;
  using span_type = typename variable_mapping<VariableList, MdSpan>::span_type;
  using extents_type = typename mdspan_type::extents_type;
  using mapping_type = typename mdspan_type::mapping;

  /// @{
  /// \name Constructors

  basic_variable_view() = default;
  basic_variable_view(const basic_variable_view&) = default;
  basic_variable_view& operator=(const basic_variable_view&) = default;
  basic_variable_view(basic_variable_view&&) = default;
  basic_variable_view& operator=(basic_variable_view&&) = default;

  template <typename OtherMdSpan,
            typename = std::enable_if_t<
                std::is_constructible<MdSpan, OtherMdSpan>::value>>
  basic_variable_view(basic_variable_view<VariableList, OtherMdSpan> other)
      : VariableList(other.get_variable_list()),
        mapping_type(other.get_extents()), m_span{other.span()} {
    assert(static_size(get_variable_list(), get_extents()) <= m_span.size());
  }

  basic_variable_view(VariableList list, span_type span,
                      extents_type extents = extents_type())
      : VariableList(list), mapping_type(extents), m_span{span} {
    assert(static_size(get_variable_list(), get_extents()) <= m_span.size());
  }
  /// @}

  /// @{
  /// \name Static Observers

  /// Returns an integer size summing up the total spanned number of
  /// elements of T.
  ///
  /// This function is equivalent to `get_span().static_size()`.
  ///
  /// \return Returns either a positive integral or `dynamic_extent` which is
  /// the compile-time known number of elements of this class.
  ///
  /// \post `static_size() > 0 || static_size() == dynamic_extent`
  ///
  /// \throws Nothing.
  static constexpr std::ptrdiff_t static_size() noexcept {
    return span_type::extent;
  }

  static constexpr std::ptrdiff_t
  static_size(VariableList vars,
              extents_type extents = extents_type()) noexcept {
    return vars.size() * mapping_type(extents).required_span_size();
  }
  /// @}

  /// \name Observers

  /// Returns the run-time known size of this object.
  ///
  /// This is equivalent to `get_span().size()`
  ///
  /// \post size() > 0
  ///
  /// \throws Nothing.
  constexpr std::ptrdiff_t size() const noexcept { return span().size(); }

  /// Returns the variable_list.
  ///
  /// \throws Nothing.
  constexpr VariableList get_variable_list() const noexcept { return *this; }

  /// Returns the mapping of indices.
  ///
  /// \throws Nothing.
  constexpr mapping_type get_mapping() const noexcept { return *this; }

  /// Returns the multi-dimensional extents.
  ///
  /// \throws Nothing.
  constexpr extents_type get_extents() const noexcept {
    return get_mapping().get_extents();
  }

  /// Returns a span of all data which is being owned be this class.
  ///
  /// \throws Nothing.
  constexpr span_type span() const noexcept { return m_span; }

  /// \name Element Access

  /// Returns a reference to data associated with any variable indexed
  /// by the VariableList of this object.
  ///
  /// \throws The implementation may throw `std::out_of_range` if the given
  /// indices are not covered by `extents()`.
  ///
  /// \note In release mode this function might not bound-check the given index
  /// and NOT throw for invalid indices.
  template <typename... IndexTypes>
  constexpr variable_ref<const basic_variable_view>
  operator()(IndexTypes... is) const {
    return variable_ref<const basic_variable_view>(*this, {is...});
  }

  /// Returns a mdspan which covers data which is associated with the
  /// given variable tag.
  ///
  /// \tparam Tag VariableList::is_valid_tag<Tag>
  ///
  /// \throws The implementation may throw `std::out_of_range` if the given
  /// tag contains an invalid index.
  template <typename Tag> constexpr mdspan_type operator[](Tag tag) const {
    return variable_mapping<VariableList, MdSpan>::make_mdspan(
        span(), get_extents(), get_variable_list(), tag);
  }

private:
  span_type m_span;
};

template <typename VariableList, typename T, typename Extents>
using variable_view =
    basic_variable_view<VariableList, basic_mdspan<T, Extents>>;

} // namespace v1
} // namespace fub

#endif