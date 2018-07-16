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

#include <boost/mp11.hpp>

namespace fub {

template <typename VariableList, typename ProtoAccessor = accessor_native<>,
          typename DefaultType = double>
struct variable_accessor {
public:
  template <typename V>
  using value_type = double; // detected_or<DefaultType, value_type_t, V>;

  template <typename V>
  using element_type =
      std::conditional_t<std::is_const<V>::value,
                         std::add_const_t<value_type<V>>, value_type<V>>;

  template <typename V>
  using accessor = typename ProtoAccessor::template rebind<element_type<V>>;

  template <typename V> using pointer = typename accessor<V>::pointer;

  template <typename V>
  using const_pointer = typename accessor<std::add_const_t<V>>::pointer;

  template <typename V> using reference = typename accessor<V>::reference;

private:
  template <template <typename...> class List, typename... Vars>
  static constexpr std::tuple<
      std::array<value_type<Vars>, variable_traits<Vars>::size()>...>
  value_storage_impl(const List<Vars...>&);

  template <typename V>
  static constexpr std::tuple<std::array<value_type<V>, variable_traits<V>::size()>>
  value_storage_impl(const V&);

  template <template <typename...> class List, typename... Vars>
  static constexpr std::tuple<pointer<Vars>...>
  pointer_storage_impl(const List<Vars...>&);

  template <typename V>
  static constexpr std::tuple<pointer<V>> pointer_storage_impl(const V&);

  template <template <typename...> class List, typename... Vars>
  static constexpr std::tuple<const_pointer<Vars>...>
  const_pointer_storage_impl(const List<Vars...>&);

  template <typename V>
  static constexpr std::tuple<const_pointer<V>>
  const_pointer_storage_impl(const V&);

public:
  using value_storage = decltype(
      variable_accessor::value_storage_impl(std::declval<VariableList>()));

  using pointer_storage = decltype(
      variable_accessor::pointer_storage_impl(std::declval<VariableList>()));

  using const_pointer_storage =
      decltype(variable_accessor::const_pointer_storage_impl(
          std::declval<VariableList>()));

  static constexpr pointer_storage to_pointer_storage(value_storage&);
  static constexpr const_pointer_storage to_pointer_storage(const value_storage&);

  template <typename V, std::ptrdiff_t... Es>
  static constexpr reference<V> access(const pointer_storage&, extents<Es...>,
                                       std::ptrdiff_t);
}; // namespace fub

template <typename VL, typename A, typename T>
constexpr typename variable_accessor<VL, A, T>::pointer_storage
variable_accessor<VL, A, T>::to_pointer_storage(value_storage& storage) {
  pointer_storage pointers;
  constexpr std::size_t N = std::tuple_size<pointer_storage>::value;
  boost::mp11::tuple_for_each(
      as_tuple_t<std::make_index_sequence<N>>(),
      [&](auto i) { std::get<i()>(pointers) = std::get<i()>(storage).data(); });
  return pointers;
}

template <typename VL, typename A, typename T>
constexpr typename variable_accessor<VL, A, T>::const_pointer_storage
variable_accessor<VL, A, T>::to_pointer_storage(const value_storage& storage) {
  pointer_storage pointers;
  constexpr std::size_t N = std::tuple_size<pointer_storage>::value;
  boost::mp11::tuple_for_each(
      as_tuple_t<std::make_index_sequence<N>>(),
      [&](auto i) { std::get<i()>(pointers) = std::get<i()>(storage).data(); });
  return pointers;
}

template <typename VL, typename A, typename T>
template <typename V, std::ptrdiff_t... Es>
constexpr typename variable_accessor<VL, A, T>::template reference<V>
variable_accessor<VL, A, T>::access(const pointer_storage& pointers,
                                    extents<Es...> e, std::ptrdiff_t n) {
  return *std::get<0>(pointers);
}

} // namespace fub

#endif