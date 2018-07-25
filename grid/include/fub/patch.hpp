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

#ifndef FUB_GRID_PATCH_HPP
#define FUB_GRID_PATCH_HPP

#include "fub/face.hpp"
#include "fub/grid/config.hpp"
#include "fub/mdspan.hpp"
#include "fub/tuple.hpp"
#include "fub/type_traits.hpp"
#include "fub/variables.hpp"

#include <range/v3/iterator_range.hpp>

#ifdef FUB_WITH_POLYMORPHIC_ALLOCATOR
#include <boost/container/pmr/polymorphic_allocator.hpp>
#endif

#include <cassert>

namespace fub {
////////////////////////////////////////////////////////////////////////////////
//                                                            Storage Descriptor
// {{{

/// @brief This storage descriptor uses an user defined allocator to allocate
/// create data on blocks.
template <typename ProtoAllocator = std::allocator<void>,
          typename LayoutPolicy = layout_left,
          typename ProtoAccessor = accessor_native<void>>
struct storage_descriptor {
  ProtoAllocator proto_allocator;
  ProtoAccessor proto_accessor;

  template <typename Variable, typename Extents, typename VariablesAccessor>
  static auto make_vector(Variable variable, Extents extents,
                          ProtoAllocator proto_alloc,
                          VariablesAccessor accessor) {
    using ValueType = typename VariableAccessor::template value_type<Variable>;
    using Traits = std::allocator_traits<ProtoAllocator>;
    typename Traits::template rebind_alloc<ValueType> allocator(proto_alloc);
    typename LayoutPolicy::template mapping<Extents> mapping(extents);
    const std::ptrdiff_t variable_size = accessor.size(variable);
    const std::ptrdiff_t span_size = mapping.required_span_size();
    return std::vector<T, Alloc>(span_size * variable_size, alloc);
  }

  template <typename Variables, typename Extents>
  auto make_storage(Variables variables, Extents extents) const {
    variable_accessor<Variables, ProtoAccessor> accessor(proto_accessor);
    return hana::unpack(variables, [=](auto... vars) {
      return hana::make_tuple(
          make_vector(vars, extents, proto_allocator, accessor)...);
    });
  }

  template <typename Variable, typename VariablesList, typename Storage,
            typename Extents>
  static auto view(Variable variable, Variables variables, Storage&& storage,
                   Extents extents) {
    variable_accessor<Variables, ProtoAccessor> variable_accessor(
        proto_accessor);
    auto span = variable_accessor.make_span(variable, storage, extents);
    using Accessor = typename ProtoAccessor::template rebind<T>;
    Accessor accessor(proto_accessor);
    return basic_mdspan<T, Extents, LayoutPolicy, Accessor>(span, extents,
                                                            accessor);
  }
};

#ifdef FUB_WITH_POLYMORPHIC_ALLOCATOR
using pmr_storage_descriptor =
    storage_descriptor<boost::container::pmr::polymorphic_allocator<void>>;
#endif

// }}}

////////////////////////////////////////////////////////////////////////////////
//                                                                 [class.patch]
// {{{
template <typename Variables, typename Extents,
#ifdef FUB_WITH_POLYMORPHIC_ALLOCATOR
          typename Descriptor = pmr_storage_descriptor>
#else
          typename Descriptor = storage_descriptor<>>
#endif
class patch {
public:
  using extents_type = Extents;
  using descriptor_type = Descriptor;
  using variables_tuple = as_tuple_t<Variables>;
  using storage_type =
      decltype(descriptor_type().make_storage(Extents(), variables_tuple()));

private:
  extents_type m_extents{};
  descriptor_type m_descriptor{};
  storage_type m_storage{};

public:
  /// @brief Default constructs storage by using default extents and default
  /// descriptors.
  patch()
      : m_storage{m_descriptor.make_storage(m_extents, variables_tuple{})} {}

  /// @brief Constructs storage with specified extents.
  explicit patch(const Extents& extents)
      : m_extents{extents}, m_storage{m_descriptor.make_storage(
                                m_extents, variables_tuple{})} {}

  /// @brief Constructs storage with specified descriptor alone.
  explicit patch(const Descriptor& descriptor)
      : m_descriptor{descriptor}, m_storage{m_descriptor.make_storage(
                                      m_extents, variables_tuple{})} {}

  /// @brief Constructs storage by specified extents and specified descriptor.
  patch(const Extents& extents, const Descriptor& descriptor)
      : m_extents{extents}, m_descriptor{descriptor},
        m_storage{m_descriptor.make_storage(m_extents, variables_tuple{})} {}

  /// @brief Returns the extents of this block.
  const Extents& extents() const noexcept { return m_extents; }

  auto size() const noexcept { return m_extents.size(); }

  /// @brief Returns the storage descriptor of this block.
  const Descriptor& descriptor() const noexcept { return m_descriptor; }

  /// Returns a view to mutable data associated to this block for a
  /// specified Variable type.
  ///
  /// Example:
  ///     struct Density {};
  ///     patch<std::tuple<Density>, extents<16, 16>> block{};
  ///     auto view = block.get<Density>();
  ///     static_assert(std::is_same_v<view, mdspan<double, 16, 16>);
  ///
  /// The example sets up a block containing a single dense variable called
  /// `Density`. We use the `get` method to obtain a view of data associated
  /// with `Density`over this patch. The standard implementation returns a \a
  /// span.
  template <typename Variable> auto get() {
    return m_descriptor.template view<Variable>(m_storage, m_extents,
                                                variables_tuple{});
  }

  /// @brief Returns a view to const data associated to this block for a
  /// specified variable type.
  template <typename Variable> auto get() const {
    return m_descriptor.template view<Variable>(m_storage, m_extents,
                                                variables_tuple{});
  }

  /// @brief Returns a view to mutable data associated with a specified
  /// variable type.
  ///
  /// Example:
  ///
  ///     struct Density { using value_type = double; };
  ///     static constexpr Density density{};
  ///     patch<std::tuple<Density>, Extents<16, 16>> block{};
  ///     auto view = block[density];
  ///     static_assert(std::is_same_v<view, span<double, Extents<16, 16>>);
  ///
  /// The example sets up a block containing a single dense variable called
  /// `Density`. We use the operator[] to obtain a view of data associated
  /// with `Density` over this patch. The standard implementation returns a \a
  /// span.
  template <typename Var> auto operator[](Var) noexcept { return get<Var>(); }

  template <typename Var> auto operator[](Var) const noexcept {
    return get<Var>();
  }

  bool operator==(const patch& other) const noexcept {
    return m_extents == other.m_extents && m_descriptor == other.m_descriptor &&
           m_storage == other.m_storage;
  }

  friend void swap(patch& b1, patch& b2) noexcept {
    using std::swap;
    swap(b1.m_extents, b2.m_extents);
    swap(b1.m_descriptor, b2.m_descriptor);
    swap(b1.m_storage, b2.m_storage);
  }

private:
  template <typename Archive>
  friend void serialize(Archive& archive, patch& p, unsigned) {
    archive& p.m_extents;
    for_each_tuple_element([&](auto& storage) { archive& storage; },
                           p.m_storage);
  }
};

template <typename Variables, typename Extents>
bool operator!=(const patch<Variables, Extents>& b1,
                const patch<Variables, Extents>& b2) noexcept {
  return !(b1 == b2);
}

template <typename Variables, typename Extents, typename Desc>
struct add_flux<patch<Variables, Extents, Desc>> {
  using type = patch<add_flux_t<Variables>, Extents, Desc>;
};

template <typename Variables, typename Extents, typename Desc>
struct add_simd<patch<Variables, Extents, Desc>> {
  using type = patch<add_simd_t<Variables>, Extents, Desc>;
};

template <typename Vars, typename Extents, typename Desc>
auto make_patch(const Vars&, const Extents& extents, const Desc& desc) {
  return patch<Vars, Extents, Desc>(extents, desc);
}

template <typename Vars, typename Extents>
auto make_patch(const Vars&, const Extents& extents) {
  return patch<Vars, Extents>(extents);
}

template <typename...> struct is_patch : std::false_type {};
template <typename T, typename E, typename D>
struct is_patch<patch<T, E, D>> : std::true_type {};

template <axis Axis, int Width, typename Vars, typename Extents, typename Desc>
auto change_extents(const patch<Vars, Extents, Desc>& p) {
  auto changed_extents = replace_extent<as_int(Axis), Width>(p.extents());
  return make_patch(Vars(), changed_extents, p.descriptor());
}

// }}}

} // namespace fub

#endif // !BLOCK_HPP
