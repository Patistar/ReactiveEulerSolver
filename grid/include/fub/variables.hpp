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

#ifndef FUB_VARIABLES_HPP
#define FUB_VARIABLES_HPP

#include "fub/tuple.hpp"
#include "fub/type_traits.hpp"

#include <meta/meta.hpp>

#include <boost/hana/size.hpp>
#include <boost/hana/tuple.hpp>

#include <string>
#include <typeinfo>
#include <vector>

namespace fub {

///////////////////////////////////////////////////////////////////////////////
//                                                            [variable_traits]
// {{{

/// Example:
// clang-format off
///     // Construct a simple variable
///     struct Density {};
///     REQUIRE(variable_traits<Density>::as_list() == meta::list<Density>());
///     REQUIRE(variable_traits<Density>::size() == 1);
///
///     // Construct a list containing one variable.
///     struct MyVariableList : basic_variable_list<Density> {};
///     REQUIRE(variable_traits<MyVariableList>::as_list() == meta::list<Density>()); 
///     REQUIRE(variable_traits<MyVariableList>::size() == 1);
// clang-format on
template <typename Variable> class variable_traits {
  /////////////////////////////////////////////////////////////////////////////
  //                                                     [variable_traits.name]
private:
  template <typename V> using name_member_function_t = decltype(V::name());

  // Returns Variable::name() because it exists.
  static const char* name_dispatch(std::true_type) noexcept {
    return Variable::name();
  }

  // Returns typeid(Variable).name() because member function `Variable::name()`
  // does NOT exsist.
  static const char* name_dispatch(std::false_type) noexcept {
    return typeid(Variable).name();
  }

public:
  /// Returns a representative string literal for this variable.
  static const char* name() {
    return name_dispatch(
        is_detected_exact<const char*, name_member_function_t, Variable>());
  }

  /////////////////////////////////////////////////////////////////////////////
  //                                                 [variable_traits.children]
private:
  template <typename V>
  using as_list_member_function_t = decltype(V::as_list());

  // Returns `Variable::as_list()` as it is detected to exist.
  static constexpr auto as_list_dispatch(std::true_type) noexcept {
    return Variable::as_list();
  }

  // Returns `meta::list<Variable>` because `Variable::as_list` does not exist.
  static constexpr boost::hana::basic_tuple<Variable>
  as_list_dispatch(std::false_type) noexcept {
    return {};
  }

public:
  /// Returns a `meta::list<Vs...>` containing all basic variables which are
  /// contained by this variable list.
  static constexpr auto as_list() noexcept {
    return as_list_dispatch(is_detected<as_list_member_function_t, Variable>());
  }

  /// Returns the number of variables contained by this variable list.
  static constexpr std::size_t size() noexcept {
    return boost::hana::size(as_list())();
  }
};

template <typename... Vs> struct basic_variable_list {
  static_assert(
      conjunction<bool_constant<variable_traits<Vs>::size() == 1>...>::value,
      "Do not merge variable lists.");
  static constexpr auto as_list() {
    return boost::hana::make_basic_tuple(boost::hana::type_c<Vs>...);
  }

  static_assert(boost::hana::size(basic_variable_list::as_list())() ==
                    sizeof...(Vs),
                "as_list does not have the correct size!");
};

// }}}

} // namespace fub

#endif
