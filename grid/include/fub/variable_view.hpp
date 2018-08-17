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
#include "fub/variable_row.hpp"

namespace fub {
inline namespace v1 {
template <typename T>
using get_extents_t = decltype(std::declval<T>().get_extents());

/// This is view on contiguous data but indexed by a variable list and a uniform
/// multi-dimensional mapping.
/// \nosubgrouping
template <typename VariableList, typename MdSpan>
class basic_variable_view : VariableList, MdSpan::mapping {
public:
  /// \name Types

  using mdspan_type = MdSpan;
  using value_type = typename mdspan_type::value_type;
  using element_type = typename mdspan_type::element_type;
  using accessor_type = typename mdspan_type::accessor;
  using span_type = typename MdSpan::span_type;
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

  template <typename Patch>
  basic_variable_view(Patch&& patch)
      : VariableList(patch.get_variable_list()),
        mapping_type(patch.get_extents()), m_span{patch.span()} {
    assert(static_size(get_variable_list(), get_extents()) <= m_span.size());
  }
  /// @}

  /// @{
  /// \name Static Observers

  static constexpr std::ptrdiff_t rank() noexcept {
    return mdspan_type::rank();
  }

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
  template <typename... IndexTypes,
            typename = std::enable_if_t<(sizeof...(IndexTypes) == rank())>,
            typename = std::enable_if_t<conjunction<
                std::is_convertible<IndexTypes, std::ptrdiff_t>...>::value>>
  constexpr variable_ref<VariableList, element_type, accessor_type>
  operator()(IndexTypes... is) const {
    std::ptrdiff_t index = get_mapping()(is...);
    return {get_variable_list(), m_span.data() + index,
            get_mapping().required_span_size(), m_span.get_accessor()};
  }

  template <typename... IndexTypes>
  constexpr variable_ref<VariableList, element_type, accessor_type>
  operator()(std::array<std::ptrdiff_t, rank()> idx) const {
    std::ptrdiff_t index = fub::apply(get_mapping(), idx);
    return {get_variable_list(), m_span.data() + index,
            get_mapping().required_span_size(), m_span.get_accessor()};
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

  constexpr variable_iterator<VariableList, element_type, accessor_type>
  begin() const noexcept {
    return {get_variable_list(), m_span.data(),
            get_mapping().required_span_size(), m_span.get_accessor()};
  }

  constexpr variable_iterator<VariableList, element_type, accessor_type>
  end() const noexcept {
    std::ptrdiff_t size = get_mapping().required_span_size();
    return {get_variable_list(), m_span.data() + size, size,
            m_span.get_accessor()};
  }

private:
  span_type m_span;
};

template <typename VariableList, typename T, typename Extents,
          typename Accessor = accessor_native<T>>
using variable_view =
    basic_variable_view<VariableList,
                        basic_mdspan<T, Extents, layout_left, Accessor>>;

template <typename VariableList, typename MdSpan>
constexpr auto rows(basic_variable_view<VariableList, MdSpan> view) noexcept {
  using Span = typename MdSpan::span_type;
  return basic_variable_row_range<VariableList, Span>{
      basic_variable_row_iterator<VariableList, Span>{
          view.get_variable_list(), view.span(),
          view.get_mapping().required_span_size(),
          row(view.get_extents()).extent(0)}};
}

template <typename L, typename M, typename R, typename FaceFeedback>
FaceFeedback for_each_face(L left, M mid, R right, FaceFeedback feedback) {
  return for_each_face(0, left, mid, right, feedback);
}

template <typename L, typename M, typename R, typename FaceFeedback>
FaceFeedback for_each_face(int dim, L left, M mid, R right,
                           FaceFeedback feedback) {
  for_each_index(mid.get_mapping(), [&](auto... is) {
    constexpr int Rank = mid.rank();
    std::array<std::ptrdiff_t, Rank> iL{{is...}};
    std::array<std::ptrdiff_t, Rank> iR = shift(iL, dim, 1);
    if (iL[dim] == 0) {
      const int last_extent_left = left.get_extents().extent(0) - 1;
      std::array<std::ptrdiff_t, Rank> iLL = replace(iL, 0, last_extent_left);
      feedback(left(iLL), mid(iL));
    }
    if (iR[dim] < mid.get_extents().extent(dim)) {
      feedback(mid(iL), mid(iR));
    } else {
      std::array<std::ptrdiff_t, Rank> iR_ = replace(iR, dim, 0);
      feedback(mid(iL), right(iR_));
    }
  });
  return feedback;
}
} // namespace v1
} // namespace fub

#endif