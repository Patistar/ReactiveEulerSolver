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

#ifndef FUB_GRID_VARIABLE_ACCESSOR_BASIC_HPP
#define FUB_GRID_VARIABLE_ACCESSOR_BASIC_HPP

#include "fub/extents.hpp"
#include "fub/quantities_ref.hpp"
#include "fub/variables.hpp"

#include <boost/hana/experimental/types.hpp>
#include <boost/hana/ext/std/integral_constant.hpp>
#include <boost/hana/for_each.hpp>
#include <boost/hana/if.hpp>
#include <boost/hana/range.hpp>
#include <boost/hana/type.hpp>

#include <array>

namespace fub {
namespace hana = boost::hana;
namespace detail {
template <typename Variable, typename ProtoAccessor, typename DefaultType>
struct variable_accessor_impl {
  template <typename V> using value_type_t = typename V::value_type;

  template <typename V,
            typename = std::enable_if_t<std::is_same<V, Variable>::value>>
  static constexpr auto get_value_type(hana::basic_type<V>) {
    return hana::type_c<detected_or_t<DefaultType, value_type_t, Variable>>;
  }

  template <typename V,
            typename = std::enable_if_t<std::is_same<V, Variable>::value>>
  static constexpr auto get_element_type(hana::basic_type<V>) {
    constexpr auto is_const = hana::trait<std::is_const>;
    constexpr auto add_const = hana::metafunction<std::add_const>;
    constexpr auto value_type = get_value_type(hana::type_c<V>);
    return hana::if_(is_const(hana::type_c<V>), add_const(value_type),
                     value_type);
  }

  template <typename V,
            typename = std::enable_if_t<std::is_same<V, Variable>::value>>
  static constexpr auto get_accessor(hana::basic_type<V>) {
    using element_type = typename decltype(get_element_type(hana::type_c<V>))::type;
    using Result = typename ProtoAccessor::template rebind<element_type>;
    return hana::type_c<Result>;
  }

  template <typename V,
            typename = std::enable_if_t<std::is_same<V, Variable>::value>>
  static constexpr auto get_value_storage(hana::basic_type<V>) {
    using value_type = typename decltype(get_value_type(hana::type_c<V>))::type;
    return hana::if_(hana::bool_c<variable_traits<V>::size() == 1>, value_type,
                     hana::type_c<std::array<value_type, size>>);
  }

  template <typename V,
            typename = std::enable_if_t<std::is_same<V, Variable>::value>>
  static constexpr auto get_pointer_storage(hana::basic_type<V>) {
    using accessor = typename decltype(get_accessor(hana::type_c<V>))::type;
    return hana::type_c<typename accessor::pointer>;
  }

  template <typename V,
            typename = std::enable_if_t<std::is_same<V, Variable>::value>>
  static constexpr auto get_const_pointer_storage(hana::basic_type<V>) {
    using accessor = typename decltype(get_accessor(hana::type_c<V>))::type;
    return hana::type_c<typename accessor::const_pointer>;
  }

  using value_storage =
      typename decltype(get_value_storage(hana::type_c<Variable>))::type;

  using pointer_storage =
      typename decltype(get_pointer_storage(hana::type_c<Variable>))::type;

  using const_pointer_storage = typename decltype(
      get_const_pointer_storage(hana::type_c<Variable>))::type;

  /// Transforms a value storage tuple into pointer storage tuple where each
  /// tuple entry points to the value tuple entry.
  /// @{
  static constexpr pointer_storage to_pointer_storage(std::true_type,
                                                      value_storage& value) {
    return std::addressof(value);
  }

  static constexpr pointer_storage to_pointer_storage(std::false_type,
                                                      value_storage& value) {
    return value.data();
  }

  static constexpr pointer_storage to_pointer_storage(value_storage& value) {
    return to_pointer_storage(bool_c<variable_traits<Variable>::size() == 1>,
                              value);
  }
  /// @}

  /// Transforms a const value storage tuple into pointer to const storage tuple
  /// where each tuple entry points to the value tuple entry.
  /// @{
  static constexpr const_pointer_storage
  to_pointer_storage(std::true_type, const value_storage& value) {
    return std::addressof(value);
  }

  static constexpr const_pointer_storage
  to_pointer_storage(std::false_type, const value_storage& value) {
    return value.data();
  }

  static constexpr const_pointer_storage
  to_pointer_storage(const value_storage&) {
    return to_pointer_storage(bool_c<variable_traits<Variable>::size() == 1>,
                              value);
  }
  /// @}

  /// Returns a reference to a variable value with specified storage, total
  /// extents and a given offset.
  template <typename V, std::ptrdiff_t... Es,
            typename = std::enable_if_t<std::is_same<V, Variable>::value>>
  static constexpr decltype(auto)
  access(hana::basic_type<V>, const pointer_storage& pointer,
         const fub::extents<Es...>& extents, std::ptrdiff_t offset) {
    const std::ptrdiff_t size = fub::size(extents);
    assert(0 <= offset && offset < size);
    typename decltype(get_accessor(hana::type_c<V>))::type accessor {}
    return accessor.access(pointer, size, offset);
  }
};
} // namespace detail

template <typename Variable, typename ProtoAccessor = accessor_native<>,
          typename DefaultType = double>
struct variable_accessor
    : private detail::variable_accessor_impl<Variable, ProtoAccessor,
                                             DefaultType> {
private:
  using base =
      detail::variable_accessor_impl<Variable, ProtoAccessor, DefaultType>;

public:
  // Variable dependent Typedefs

  template <typename V>
  using value_type =
      typename decltype(base::get_value_type(hana::type_c<V>))::type;

  template <typename V>
  using element_type =
      typename decltype(base::get_element_type(hana::type_c<V>))::type;

  template <typename V>
  using accessor =
      typename decltype(base::get_accessor_type(hana::type_c<V>))::type;

  template <typename V> using pointer = typename accessor<V>::pointer;

  template <typename V>
  using const_pointer = typename accessor<std::add_const_t<V>>::pointer;

  template <typename V> using reference = typename accessor<V>::reference;

  template <typename V>
  using const_reference = typename accessor<std::add_const_t<V>>::reference;

  // Storage related Typedefs

  using value_storage = typename decltype(
      base::get_value_storage(std::declval<VariableList>()))::type;

  using pointer_storage = typename decltype(
      base::get_pointer_storage(std::declval<VariableList>()))::type;

  using const_pointer_storage = typename decltype(
      base::const_pointer_storage(std::declval<VariableList>()))::type;

  /// Transforms a value storage tuple into pointer storage tuple where each
  /// tuple entry points to the value tuple entry.
  static constexpr pointer_storage to_pointer_storage(value_storage&);

  /// Transforms a const value storage tuple into pointer to const storage tuple
  /// where each tuple entry points to the value tuple entry.
  static constexpr const_pointer_storage
  to_pointer_storage(const value_storage&);

  /// Returns a reference to a variable value with specified storage, total
  /// extents and a given offset.
  template <typename V, std::ptrdiff_t... Es>
  static constexpr reference<V>
  access(hana::basic_tuple<V> variable, const pointer_storage& storage,
         const extents<Es...>& extents, std::ptrdiff_t offset);
}; // namespace fub

template <typename VL, typename A, typename T>
constexpr typename variable_accessor<VL, A, T>::pointer_storage
variable_accessor<VL, A, T>::to_pointer_storage(value_storage& storage) {
  pointer_storage pointers;
  constexpr std::size_t N = std::tuple_size<pointer_storage>::value;
  hana::for_each(hana::make_range(hana::int_c<0>, hana::int_c<N>), [&](auto i) {
    std::get<i()>(pointers) = std::get<i()>(storage).data();
  });
  return pointers;
}

template <typename VL, typename A, typename T>
constexpr typename variable_accessor<VL, A, T>::const_pointer_storage
variable_accessor<VL, A, T>::to_pointer_storage(const value_storage& storage) {
  pointer_storage pointers;
  constexpr std::size_t N = std::tuple_size<pointer_storage>::value;
  hana::for_each(hana::make_range(hana::int_c<0>, hana::int_c<N>), [&](auto i) {
    std::get<i()>(pointers) = std::get<i()>(storage).data();
  });
  return pointers;
}

template <typename VL, typename A, typename T>
template <typename V, std::ptrdiff_t... Es>
constexpr typename variable_accessor<VL, A, T>::template reference<V>
variable_accessor<VL, A, T>::access(hana::type<V>,
                                    const pointer_storage& pointers,
                                    const extents<Es...>& total_extents,
                                    std::ptrdiff_t offset) {
  constexpr auto index =
      hana::index_if(hana_types(), hana::equal.to(hana::type_c<V>));
  static_assert(index != hana::nothing,
                "Could not find specified Variable in the VariableList.");
  return *std::get<*index>(pointers);
}

} // namespace fub

#endif