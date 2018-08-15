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

#ifndef FUB_GRID_VARIABLE_DATA_HPP
#define FUB_GRID_VARIABLE_DATA_HPP

#include "fub/simd.hpp"
#include "fub/variable_mapping.hpp"
#include "fub/variable_row.hpp"

#include <boost/align/aligned_allocator.hpp>

#include <array>

namespace fub {
inline namespace v1 {
namespace detail {
template <std::ptrdiff_t N, typename VariableList, typename MdSpan,
          typename Allocator>
struct basic_variable_data_storage : VariableList, MdSpan::mapping {
  using value_type = typename Allocator::value_type;
  using mapping_type = typename MdSpan::mapping;

  basic_variable_data_storage() = default;

  basic_variable_data_storage(VariableList list, mapping_type mapping,
                              Allocator alloc = Allocator())
      : VariableList(list), mapping_type(mapping) {}

  constexpr fub::span<value_type, N> span() noexcept {
    return make_span(m_data);
  }

  constexpr fub::span<const value_type, N> span() const noexcept {
    return make_span(m_data);
  }

  constexpr const mapping_type& get_mapping() const noexcept { return *this; }

  constexpr const VariableList& get_variable_list() const noexcept {
    return *this;
  }

  static constexpr int alignment = Vc::Vector<value_type>::MemoryAlignment;

  alignas(alignment) std::array<value_type, N> m_data;
};

template <typename VariableList, typename MdSpan, typename Allocator>
struct basic_variable_data_storage<dynamic_extent, VariableList, MdSpan,
                                   Allocator> : VariableList,
                                                MdSpan::mapping {
  using value_type = typename Allocator::value_type;
  using mapping_type = typename MdSpan::mapping;

  basic_variable_data_storage() = default;

  basic_variable_data_storage(VariableList list, mapping_type mapping,
                              Allocator alloc = Allocator())
      : VariableList(list), mapping_type(mapping),
        m_data(list.size() * mapping.required_span_size(), alloc) {
    assert(list.size() > 0);
    assert(mapping.required_span_size() > 0);
  }

  fub::span<value_type> span() noexcept {
    return {m_data.data(), static_cast<std::ptrdiff_t>(m_data.size())};
  }

  fub::span<const value_type> span() const noexcept {
    return {m_data.data(), static_cast<std::ptrdiff_t>(m_data.size())};
  }

  const mapping_type& get_mapping() const noexcept { return *this; }

  const VariableList& get_variable_list() const noexcept { return *this; }

  std::vector<value_type, Allocator> m_data;
};
} // namespace detail

/// This class provides data storage for a given variable list and uniform
/// specified extents.
///
/// If the total size of data is a compile-time known this class will never use
/// the Allocator to allocate memory but will use the automatic storage instead.
///
/// \note Note that the MdSpan contains information about the extents, layout
/// mapping and accessor policy for each variable.
template <typename VariableList, typename MdSpan,
          typename Allocator = boost::alignment::aligned_allocator<
              typename MdSpan::value_type,
              Vc::Vector<typename MdSpan::value_type>::MemoryAlignment>>
class basic_variable_data {
  using variable_map = variable_mapping<VariableList, MdSpan>;

public:
  /// \name Static Observers

  /// Returns either a positive integral or `dynamic_extent` which is
  /// the compile-time known number of elements of this class.
  ///
  /// This function is equivalent to `get_span().static_size()`.
  ///
  /// \throws Nothing.
  static constexpr std::ptrdiff_t static_size() noexcept {
    return variable_map::static_size();
  }

  /// \name Types

  using value_type = typename Allocator::value_type;
  using allocator_type = Allocator;
  using accessor_type = typename variable_map::accessor_type;
  using span_type = typename variable_map::span_type;
  using const_span_type = typename variable_map::const_span_type;
  using mdspan_type = typename variable_map::mdspan_type;
  using const_mdspan_type = typename variable_map::const_mdspan_type;
  using extents_type = typename variable_map::extents_type;
  using mapping_type = typename variable_map::mapping_type;
  using storage_type =
      detail::basic_variable_data_storage<basic_variable_data::static_size(),
                                          VariableList, MdSpan, Allocator>;

  /// \name Constructors

  basic_variable_data() = default;
  basic_variable_data(const basic_variable_data&) = default;
  basic_variable_data& operator=(const basic_variable_data&) = default;
  basic_variable_data(basic_variable_data&&) = default;
  basic_variable_data& operator=(basic_variable_data&&) = default;

  explicit basic_variable_data(VariableList list, extents_type extents = extents_type(),
                      allocator_type alloc = allocator_type())
      : m_storage(list, extents, alloc) {}

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
  constexpr VariableList get_variable_list() const noexcept {
    return m_storage.get_variable_list();
  }

  /// Returns the mapping of indices.
  ///
  /// \throws Nothing.
  constexpr mapping_type get_mapping() const noexcept {
    return m_storage.get_mapping();
  }

  /// Returns the multi-dimensional extents.
  ///
  /// \throws Nothing.
  constexpr extents_type get_extents() const noexcept {
    return get_mapping().get_extents();
  }

  /// Returns a span of all data which is being owned be this class.
  ///
  /// \throws Nothing.
  /// @{
  constexpr span_type span() noexcept { return m_storage.span(); }
  constexpr const_span_type span() const noexcept { return m_storage.span(); }
  /// @}

  /// \name Element Access

  /// Returns a reference to data associated with any variable indexed
  /// by the VariableList of this object.
  ///
  /// \throws The implementation may throw `std::out_of_range` if the given
  /// indices are not covered by `extents()`.
  ///
  /// \note In release mode this function might not check given index range and
  /// NOT throw on invalid indices.
  /// @{
  template <typename... IndexTypes>
  constexpr variable_ref<VariableList, value_type>
  operator()(IndexTypes... is) {
    std::ptrdiff_t index = get_mapping()(is...);
    return {get_variable_list(), m_storage.span().data() + index,
            get_mapping().required_span_size()};
  }

  template <typename... IndexTypes>
  constexpr variable_ref<VariableList, const value_type>
  operator()(IndexTypes... is) const {
    std::ptrdiff_t index = get_mapping()(is...);
    return {get_variable_list(), m_storage.data() + index,
            get_mapping().required_span_size()};
  }
  /// @}

  /// Returns a mdspan which covers data which is associated with the
  /// given variable tag.
  ///
  /// \requires VariableList::is_valid_tag<Tag>
  ///
  /// \throws Nothing.
  /// @{
  template <typename Tag> constexpr mdspan_type operator[](Tag tag) {
    return variable_map::make_mdspan(span(), get_extents(), get_variable_list(),
                                     tag);
  }
  template <typename Tag>
  constexpr const_mdspan_type operator[](Tag tag) const {
    return variable_map::make_mdspan(span(), get_extents(), get_variable_list(),
                                     tag);
  }
  /// @}

  variable_iterator<VariableList, value_type, accessor_type> begin() noexcept {
    return {get_variable_list(), span().data(),
            get_mapping().required_span_size(), span().get_accessor()};
  }

  variable_iterator<VariableList, value_type, accessor_type> end() noexcept {
    std::ptrdiff_t size = get_mapping().required_span_size();
    return {get_variable_list(), span().data() + size, size,
            span().get_accessor()};
  }

  variable_iterator<VariableList, const value_type, accessor_type> begin() const
      noexcept {
    return {get_variable_list(), span().data(),
            get_mapping().required_span_size(), span().get_accessor()};
  }

  variable_iterator<VariableList, const value_type, accessor_type> end() const
      noexcept {
    std::ptrdiff_t size = get_mapping().required_span_size();
    return {get_variable_list(), span().data() + size, size,
            span().get_accessor()};
  }

private:
  storage_type m_storage;
};

template <typename VariableList, typename T, int Rank>
using variable_data =
    basic_variable_data<VariableList, basic_mdspan<T, dynamic_extents_t<Rank>>,
                        std::allocator<T>>;

template <typename VariableList, typename MdSpan, typename Allocator>
constexpr auto
rows(basic_variable_data<VariableList, MdSpan, Allocator>& view) noexcept {
  using Span = typename MdSpan::span_type;
  return basic_variable_row_range<VariableList, Span>{
      basic_variable_row_iterator<VariableList, Span>{
          view.get_variable_list(), view.span(),
          view.get_mapping().required_span_size(),
          row(view.get_extents()).extent(0)},
      rows(view.get_extents())};
}

template <typename VariableList, typename MdSpan, typename Allocator>
constexpr auto rows(
    const basic_variable_data<VariableList, MdSpan, Allocator>& view) noexcept {
  using Span = typename basic_variable_data<VariableList, MdSpan,
                                            Allocator>::const_span_type;
  return basic_variable_row_range<VariableList, Span>{
      basic_variable_row_iterator<VariableList, Span>{
          view.get_variable_list(), view.span(),
          view.get_mapping().required_span_size(),
          row(view.get_extents()).extent(0)},
      rows(view.get_extents())};
}

} // namespace v1
} // namespace fub

#endif