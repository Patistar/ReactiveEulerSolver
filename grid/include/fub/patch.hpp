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
template <typename D, typename E>
using has_make_storage_member_function_t =
    decltype(std::declval<D&>().make_storage(std::declval<const E&>()));

template <typename D, typename E>
struct has_make_storage_member_function
    : is_detected<has_make_storage_member_function_t, D, E> {};

template <typename D, typename S, typename E>
using has_view_member_function_t = decltype(
    std::declval<D&>().view(std::declval<S>(), std::declval<const E&>()));

template <typename D, typename S, typename E>
struct has_view_member_function
    : is_detected<has_view_member_function_t, D, S, E> {};

template <typename Descriptor> struct descriptor_traits {
  using descriptor_type = Descriptor;

  template <typename Extents, typename... Vars>
  static auto make_storage(Descriptor& descriptor, const Extents& extents,
                           const std::tuple<Vars...>& variables) {
    return descriptor.make_storage(extents, variables);
  }

  template <typename Extents, typename Vars>
  using storage_type = decltype(
      make_storage(std::declval<Descriptor&>(), std::declval<const Extents&>(),
                   std::declval<const as_tuple_t<Vars>&>()));

  template <typename Variable, typename Extents, typename Vars>
  static auto view(const Descriptor& /* descriptor */,
                   storage_type<Extents, Vars>& storage, const Extents& extents,
                   const Vars&) {
    using pointer = typename variable_traits<Variable>::pointer;
    using index = variable_find_index<Variable, Vars>;
    pointer data = std::get<index::value>(storage).data();
    return mdspan<typename variable_traits<Variable>::value_type, Extents>{
        data, extents};
  }

  template <typename Variable, typename Extents, typename Vars>
  static auto view(const Descriptor& /* descriptor */,
                   const storage_type<Extents, Vars>& storage,
                   const Extents& extents, const Vars&) {
    using index = variable_find_index<Variable, Vars>;
    auto data = std::get<index::value>(storage).data();
    return mdspan<const typename variable_traits<Variable>::value_type,
                  Extents>{data, extents};
  }

  template <typename Variable, typename Extents, typename Variables>
  using view_type = decltype(
      view(std::declval<const Variable&>(),
           std::declval<storage_type<Extents, Variables>&>(),
           std::declval<const Extents&>(), std::declval<const Variables&>()));

  template <typename Variable, typename Extents, typename Variables>
  using const_view_type = decltype(
      view(std::declval<const Variable&>(),
           std::declval<const storage_type<Extents, Variables>&>(),
           std::declval<const Extents&>(), std::declval<const Variables&>()));
};

/// @brief This storage descriptor tries to use an `std::array<T, N>` whenever
/// all extents are compile-time known.
struct automatic_storage_descriptor {
  template <typename Extents, typename... Vars>
  static auto make_storage(const Extents& e, const std::tuple<Vars...>& vars,
                           std::true_type) noexcept {
    return fub::apply(
        [](auto... vs) {
          constexpr index Size = Extents().size();
          auto make_array = [](auto v) {
            // We do not initialize the array here. Just create the storage.
            std::array<typename variable_traits<decltype(v)>::value_type, Size>
                array;
            return array;
          };
          return std::make_tuple(make_array(vs)...);
        },
        vars);
  }

  template <typename Extents, typename... Vars>
  static auto make_storage(const Extents& extents,
                           const std::tuple<Vars...>& vars, std::false_type) {
    return fub::apply(
        [&](auto... vs) {
          auto make_vector = [](auto v, const Extents& extents) {
            const auto size = variable_traits<decltype(v)>::size(extents);
            return std::vector<
                typename variable_traits<decltype(v)>::value_type>(size);
          };
          return std::make_tuple(make_vector(vs, extents)...);
        },
        vars);
  }

  template <typename Extents, typename... Vars>
  static auto make_storage(
      const Extents& extents,
      const std::tuple<Vars...>& vars) noexcept(Extents::rank_dynamic == 0) {
    return make_storage(extents, vars,
                        bool_c < Extents::rank_dynamic == 0 &&
                            (sizeof...(Vars) < 20) >);
  }
};

/// @brief This storage descriptor uses an user defined allocator to allocate
/// create data on blocks.
template <typename Allocator = std::allocator<void>>
struct dynamic_storage_descriptor {
  using proto_allocator = Allocator;
  proto_allocator allocator;

  template <typename Extents, typename... Vars>
  auto make_storage(const Extents& extents,
                    const std::tuple<Vars...>& vars) const {
    return fub::apply(
        [&extents, alloc = allocator](auto... vs) {
          auto make_pmr_vector = [](auto v, const Extents& extents,
                                    const proto_allocator& alloc) {
            using Var = decltype(v);
            using T = typename variable_traits<Var>::value_type;
            using Traits = std::allocator_traits<proto_allocator>;
            using Alloc = typename Traits::template rebind_alloc<T>;
            const std::ptrdiff_t size = variable_traits<Var>::size(extents);
            return std::vector<T, Alloc>(size, Alloc{alloc});
          };
          return std::make_tuple(make_pmr_vector(vs, extents, alloc)...);
        },
        vars);
  }
};

#ifdef FUB_WITH_POLYMORPHIC_ALLOCATOR
using pmr_storage_descriptor = dynamic_storage_descriptor<
    boost::container::pmr::polymorphic_allocator<void>>;
#endif

// }}}

////////////////////////////////////////////////////////////////////////////////
//                                                                 [class.patch]
// {{{
template <typename Variables, typename Extents,
          typename Descriptor = pmr_storage_descriptor>
class patch {
  using Traits = descriptor_traits<Descriptor>;

public:
  using extents_type = Extents;
  using descriptor_type = Descriptor;
  using variables_tuple = as_tuple_t<Variables>;
  using storage_type =
      typename Traits::template storage_type<Extents, Variables>;

private:
  extents_type m_extents{};
  descriptor_type m_descriptor{};
  storage_type m_storage{};

public:
  /// @brief Default constructs storage by using default extents and default
  /// descriptors.
  patch()
      : m_storage{
            Traits::make_storage(m_descriptor, m_extents, variables_tuple{})} {}

  /// @brief Constructs storage with specified extents.
  explicit patch(const Extents& extents)
      : m_extents{extents}, m_storage{Traits::make_storage(
                                m_descriptor, m_extents, variables_tuple{})} {}

  /// @brief Constructs storage with specified descriptor alone.
  explicit patch(const Descriptor& descriptor)
      : m_descriptor{descriptor}, m_storage{Traits::make_storage(
                                      m_descriptor, m_extents,
                                      variables_tuple{})} {}

  /// @brief Constructs storage by specified extents and specified descriptor.
  patch(const Extents& extents, const Descriptor& descriptor)
      : m_extents{extents}, m_descriptor{descriptor},
        m_storage{
            Traits::make_storage(m_descriptor, m_extents, variables_tuple{})} {}

  /// @brief Returns the extents of this block.
  const Extents& extents() const noexcept { return m_extents; }

  auto size() const noexcept { return m_extents.size(); }

  /// @brief Returns the storage descriptor of this block.
  const Descriptor& descriptor() const noexcept { return m_descriptor; }

  /// @brief Returns a view to mutable data associated to this block for a
  /// specified Variable type.
  ///
  /// Example:
  ///
  ///     struct Density { using value_type = double; };
  ///     patch<std::tuple<Density>, extents<16, 16>> block{};
  ///     auto view = block.get<Density>();
  ///     static_assert(std::is_same_v<view, span<double, extents<16, 16>>);
  ///
  /// The example sets up a block containing a single dense variable called
  /// `Density`. We use the `get` method to obtain a view of data associated
  /// with `Density`over this patch. The standard implementation returns a \a
  /// span.
  template <typename Variable> auto get() {
    return Traits::template view<Variable>(m_descriptor, m_storage, m_extents,
                                           variables_tuple{});
  }

  /// @brief Returns a view to const data associated to this block for a
  /// specified variable type.
  template <typename Variable> auto get() const {
    return Traits::template view<Variable>(m_descriptor, m_storage, m_extents,
                                           variables_tuple{});
  }

  /// @brief Returns a view to mutable data associated with a specified
  /// variable type.
  ///
  /// Example:
  ///
  ///     struct Density { using value_type = double; };
  ///     static constexpr Density density{};
  ///     patch<Variables<Density>, Extents<16, 16>> block{};
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
    std::swap(b1.m_extents, b2.m_extents);
    std::swap(b1.m_descriptor, b2.m_descriptor);
    std::swap(b1.m_storage, b2.m_storage);
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
