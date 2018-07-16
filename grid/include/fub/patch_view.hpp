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

/// This a view type for `patch<Vs, E, S>` which provides mdspans for multiple
/// variables with uniform extents, layout and accessor policies.
///
/// This type essentially erases the StorageDescriptor from the `patch` class
/// and only points to the data located on a `patch`.
///
/// Example:
///     struct MyVariable {} variable;
///
///     // Create a fixed-size patch containing data for variable MyVariable
///     patch<meta::list<MyVariable>, extents<32, 32>> p;
///
///     // Create a view on the patch p.
///     patch_view<meta::list<MyVariable>, extents<32, 32>> view = p;
///
///     // Access the data
///     view[variable](10, 10) = 42.;
///     REQUIRE(view(10, 10)[variable] == 42.);
template <typename VariableList, typename Extents, typename Layout,
          typename AccessorPolicy>
class basic_patch_view : private Layout::template mapping<Extents> {
  static_assert(is_extents_v<Extents>,
                "Extents must be a valid extents<StaticExtents...> object.");

public:
  using mapping = typename Layout::template mapping<Extents>;
  using extents_type = typename mapping::extents_type;
  using index_type = typename extents_type::index_type;
  using variable_accessor =
      typename AccessorPolicy::template accessor<VariableList>;
  using reference = basic_quantities_ref<VariableList, AccessorPolicy>;
  using pointer_storage = typename variable_accessor::pointer_storage;

  template <typename Variable>
  using element_type = typename variable_accessor::element_type;

  template <typename Varibale>
  using mdspan =
      basic_mdspan<element_type<Variable>, Extents, layout_left, accessor>;

  basic_patch_view() = default;
  basic_patch_view(const basic_patch_view&) = default;
  basic_patch_view(basic_patch_view&&) = default;
  basic_patch_view& operator=(const basic_patch_view&) = default;
  basic_patch_view& operator=(basic_patch_view&&) = default;

  constexpr basic_patch_view(const pointer_storage& p, const Extents& e);

  template <typename Patch,
            typename std::enable_if_t<is_patch_v<
                remove_cvref_t<Patch>, VariableList, Extents>>* = nullptr>
  basic_patch_view(Patch&& patch);

  constexpr const mapping& get_mapping() const noexcept;

  constexpr const extents_type& get_extents() const noexcept;

  constexpr const variable_accessor& get_variable_accessor() const noexcept;

  constexpr std::ptrdiff_t size() const noexcept;

  template <typename Variable>
  constexpr mdspan<Variable> operator[](Variable var) const noexcept;

  template <
      typename... Index,
      typename std::enable_if_t<is_invocable_v<mapping, Index...>* = nullptr>>
  constexpr reference operator()(Index... is) const noexcept;

private:
  pointer m_pointer{};
};

template <typename VariableList, typename Extents>
using patch_view = basic_patch_view<VariableList, Extents, layout_left,
                                    variable_accessor<accessor_basic>>;

} // namespace fub

#endif // !BLOCKVIEW_HPP
