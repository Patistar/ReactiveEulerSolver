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

#ifndef FUB_QUANTITIES_HPP
#define FUB_QUANTITIES_HPP

#include "fub/variable_accessor_basic.hpp"

namespace fub {

template <typename VariableList, typename Accessor> class basic_quantities {
public:
  using accessor = Accessor;

  using variable_list = VariableList;

  template <typename V>
  using value_type = typename accessor::template value_type<V>;

  template <typename V>
  using reference = typename accessor::template reference<V>;

  template <typename V>
  using const_reference = typename accessor::template reference<std::add_const_t<V>>;

  /// Constructs a quantities tuple from specified values.
  explicit basic_quantities(const std::tuple<value_type<V>...>&, const accessor&);

  template <std::size_t I>
  decltype(auto) get() noexcept;

  template <std::size_t I>
  decltype(auto) get() const noexcept;

  template <typename V>
  decltype(auto) get() noexcept;

  template <typename V>
  decltype(auto) get() const noexcept;

  template <typename V>
  reference<V> operator[](const V&) noexcept;

  template <typename V>
  const_reference<V> operator[](const V&) const noexcept;

  static std::size_t size() const noexcept;

private:
  typename accessor::value_storage_type m_storage;
};

template <typename VariableList> 
using quantities = basic_quantities<VariableList, variable_accesor_basic>;

} // namespace fub

#endif