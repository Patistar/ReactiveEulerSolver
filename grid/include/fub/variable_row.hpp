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

#ifndef FUB_GRID_VARIABLE_ROW_HPP
#define FUB_GRID_VARIABLE_ROW_HPP

#include "fub/span.hpp"

namespace fub {
inline namespace v1 {

template <typename VariableList, typename Span>
class basic_variable_row : VariableList {
public:
  /// \name Types

  using value_type = typename Span::value_type;
  using element_type = typename Span::element_type;
  using accessor_type = typename Span::accessor;
  using span_type = Span;

  /// @{
  /// \name Constructors

  basic_variable_row() = default;
  basic_variable_row(const basic_variable_row&) = default;
  basic_variable_row& operator=(const basic_variable_row&) = default;
  basic_variable_row(basic_variable_row&&) = default;
  basic_variable_row& operator=(basic_variable_row&&) = default;

  template <typename OtherSpan,
            typename =
                std::enable_if_t<std::is_constructible<Span, OtherSpan>::value>>
  basic_variable_row(basic_variable_row<VariableList, OtherSpan> other)
      : VariableList(other.get_variable_list()), m_span{other.span()},
        m_stride{other.stride()}, m_length{other.length()} {}

  basic_variable_row(VariableList list, span_type span, std::ptrdiff_t stride,
                     std::ptrdiff_t length)
      : VariableList(list), m_span{span}, m_stride{stride}, m_length{length} {}

  /// @}

  /// @{
  /// \name Static Observers

  /// Returns an integer size summing up the total spanned number of
  /// elements of T.
  ///
  /// This function is equivalent to `span().static_size()`.
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
  /// @}

  /// \name Observers

  /// Returns the run-time known size of this object.
  ///
  /// This is equivalent to `span().size()`
  ///
  /// \post size() > 0
  ///
  /// \throws Nothing.
  constexpr std::ptrdiff_t size() const noexcept { return span().size(); }
  constexpr std::ptrdiff_t stride() const noexcept { return m_stride; }
  constexpr std::ptrdiff_t length() const noexcept { return m_length; }

  /// Returns the variable_list.
  ///
  /// \throws Nothing.
  constexpr VariableList get_variable_list() const noexcept { return *this; }

  /// Returns a span of all data which is being owned be this class.
  ///
  /// \throws Nothing.
  constexpr span_type span() const noexcept { return m_span; }

  constexpr auto next() {
    auto next_span = drop(m_span, m_length);
    using next_span_t = decltype(next_span);
    return basic_variable_row<VariableList, next_span_t>{
        get_variable_list(), next_span, m_stride, m_length};
  }

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
  constexpr variable_ref<VariableList, element_type, accessor_type>
  operator()(std::ptrdiff_t index) const {
    return {get_variable_list(), m_span.data() + index, stride(),
            m_span.get_accessor()};
  }

  constexpr variable_iterator<VariableList, element_type, accessor_type>
  begin() const noexcept {
    return {get_variable_list(), m_span.data(), stride(),
            m_span.get_accessor()};
  }

  constexpr variable_iterator<VariableList, element_type, accessor_type>
  end() const noexcept {
    return {get_variable_list(), m_span.data() + length(), stride(),
            m_span.get_accessor()};
  }

private:
  span_type m_span;
  std::ptrdiff_t m_stride;
  std::ptrdiff_t m_length;
};

template <typename VariableList, typename Span>
class basic_variable_row_iterator {
public:
  /// \name Types

  using value_type = typename Span::value_type;
  using element_type = typename Span::element_type;
  using accessor_type = typename Span::accessor;
  using span_type = Span;

  /// @{
  /// \name Constructors

  basic_variable_row_iterator() = default;

  basic_variable_row_iterator(const basic_variable_row_iterator&) = default;

  basic_variable_row_iterator(basic_variable_row_iterator&&) = default;

  basic_variable_row_iterator&
  operator=(const basic_variable_row_iterator&) = default;

  basic_variable_row_iterator&
  operator=(basic_variable_row_iterator&&) = default;

  /// Constructs a row iterator that points to a row in a specific span.
  basic_variable_row_iterator(VariableList list, span_type span,
                              std::ptrdiff_t stride, std::ptrdiff_t length)
      : m_row(list, span, stride, length) {}

  /// @}

  basic_variable_row<VariableList, Span> read() const noexcept { return m_row; }

  void next() { m_row = m_row.next(); }

  bool equal_to(basic_variable_row_iterator other) const noexcept {
    return m_row.span() == other.m_row.span();
  }

private:
  basic_variable_row<VariableList, Span> m_row{};
};

template <typename VariableList, typename Span>
bool operator==(basic_variable_row_iterator<VariableList, Span> lhs,
                basic_variable_row_iterator<VariableList, Span> rhs) {
  return lhs.equal_to(rhs);
}

template <typename VariableList, typename Span>
bool operator!=(basic_variable_row_iterator<VariableList, Span> lhs,
                basic_variable_row_iterator<VariableList, Span> rhs) {
  return !(lhs == rhs);
}

template <typename VariableList, typename Span>
basic_variable_row_iterator<VariableList, Span>&
operator++(basic_variable_row_iterator<VariableList, Span>& iter) {
  iter.next();
  return iter;
}

template <typename VariableList, typename Span>
basic_variable_row_iterator<VariableList, Span>
operator++(basic_variable_row_iterator<VariableList, Span>& iter, int) {
  auto copy = iter;
  iter.next();
  return copy;
}

template <typename VariableList, typename Span>
basic_variable_row<VariableList, Span>
operator*(basic_variable_row_iterator<VariableList, Span> iter) {
  return iter.read();
}

template <typename VariableList, typename Span>
class basic_variable_row_range {
public:
  using iterator = basic_variable_row_iterator<VariableList, Span>;

  basic_variable_row_range() = default;

  basic_variable_row_range(iterator first, std::ptrdiff_t size)
      : m_first{first}, m_last{fub::next(first, size)}, m_size{size} {}

  std::ptrdiff_t size() const noexcept { return m_size; }

  iterator begin() const noexcept { return m_first; }

  iterator end() const noexcept { return m_last; }

private:
  iterator m_first{};
  iterator m_last{};
  std::ptrdiff_t m_size;
};



} // namespace v1
} // namespace fub

#endif