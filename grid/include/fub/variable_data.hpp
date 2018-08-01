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

#include "fub/variable_mapping.hpp"
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

  std::array<value_type, N> m_data;
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
template <typename VariableList, typename ValueType, typename MdSpan,
          typename Allocator = std::allocator<ValueType>>
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

  using value_type = ValueType;
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

  basic_variable_data(VariableList list, extents_type extents,
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
  constexpr variable_ref<basic_variable_data> operator()(IndexTypes... is) {
    return variable_ref<basic_variable_data>(*this, {is...});
  }
  template <typename... IndexTypes>
  constexpr variable_ref<const basic_variable_data>
  operator()(IndexTypes... is) const {
    return variable_ref<const basic_variable_data>(*this, {is...});
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

private:
  storage_type m_storage;
};

template <typename VariableList, typename T, typename Extents>
using variable_data =
    basic_variable_data<VariableList, T, basic_mdspan<T, Extents>,
                        std::allocator<T>>;

} // namespace v1
} // namespace fub

#endif