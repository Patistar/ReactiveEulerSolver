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

#include "fub/type_traits.hpp"

#include "fub/hana.hpp"

#include <boost/hana/any_of.hpp>
#include <boost/hana/equal.hpp>
#include <boost/hana/ext/std/integral_constant.hpp>
#include <boost/hana/fold.hpp>
#include <boost/hana/not_equal.hpp>
#include <boost/hana/size.hpp>
#include <boost/hana/tuple.hpp>
#include <boost/hana/type.hpp>

#include <typeinfo>
#include <vector>

namespace fub {

///////////////////////////////////////////////////////////////////////////////
//                                                            [variable_traits]
// {{{

template <typename Collection> struct collection_traits {
private:
  /// This is a helper typedef to detect if `v.as_list()` is a valid expression
  /// for a object v of type V.
  template <typename C>
  using variables_member_function_t = decltype(std::declval<C>().variables());

  // Returns `Variable::as_list()` as it is detected to exist.
  static constexpr auto variables_dispatch(std::true_type,
                                           Collection collection) noexcept {
    return collection.variables();
  }

  /// Returns `hana::tuple<Collection>` because `Collection::variables` does not
  /// exist.
  static constexpr hana::tuple<Collection>
  variables_dispatch(std::false_type, Collection collection) noexcept {
    return {collection};
  }

public:
  /// Returns a `hana::tuple<Vs...>` containing all basic variables which are
  /// contained by this collection.
  static constexpr auto variables(Collection collection) noexcept {
    return variables_dispatch(
        is_detected<variables_member_function_t, Collection>{}, collection);
  }

private:
  /// This is a helper typedef to detect if `v.as_list()` is a valid expression
  /// for a object v of type V.
  template <typename C>
  using size_member_function_t = decltype(std::declval<C>().size());

  // Returns `Variable::as_list()` as it is detected to exist.
  static constexpr auto size_dispatch(std::true_type,
                                      Collection collection) noexcept {
    return collection.size();
  }

  /// Returns `hana::tuple<Collection>` because `Collection::size` does not
  /// exist.
  static constexpr std::ptrdiff_t
  size_dispatch(std::false_type, Collection collection) noexcept {
    return hana::size(variables(collection))();
  }

public:
  /// Returns the number of variables contained by this variable list.
  static constexpr std::ptrdiff_t size(Collection collection) noexcept {
    return size_dispatch(is_detected<size_member_function_t, Collection>{},
                         collection);
  }
};

template <typename> struct is_collection : bool_constant<false> {};

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
  template <typename V>
  using name_member_function_t = decltype(std::declval<V>().name());

  // Returns Variable::name() because it exists.
  static const char* name_dispatch(std::true_type, Variable variable) noexcept {
    return Variable::name(variable);
  }

  // Returns typeid(Variable).name() because member function `Variable::name()`
  // does NOT exsist.
  static const char* name_dispatch(std::false_type, Variable) noexcept {
    return typeid(Variable).name();
  }

public:
  /// Returns a representative string literal for this variable.
  static const char* name(Variable variable) noexcept {
    return name_dispatch(
        is_detected_exact<const char*, name_member_function_t, Variable>{},
        variable);
  }

private:
  template <typename C, typename V>
  using is_iterable_in_t = decltype(
      std::declval<C>().begin() <= std::declval<C>().next(std::declval<V>()) &&
      std::declval<C>().next(std::declval<V>()) < std::declval<C>().end());

  template <typename Collection>
  static constexpr bool is_iterable_in_dispatch(std::false_type, Collection,
                                                Variable) noexcept {
    return false;
  }

  template <typename Collection>
  static constexpr bool is_iterable_in_dispatch(std::true_type,
                                                Collection collection,
                                                Variable variable) noexcept {
    return collection.begin() <= variable && variable < collection.end();
  }

public:
  template <typename Collection>
  static constexpr bool is_iterable_in(Collection collection,
                                       Variable variable) noexcept {
    return is_iterable_in_dispatch(
        is_detected<is_iterable_in_t, Collection, Variable>{}, collection,
        variable);
  }
};

template <typename... Vs>
struct static_collection : private hana::tuple<Vs...> {
  using Base = hana::tuple<Vs...>;
  using Base::Base;

  static_assert(conjunction<bool_constant<collection_traits<Vs>::size(Vs{}) <=
                                          1>...>::value,
                "Do not merge non-trivial variables.");

  static_assert(!hana::any_of(hana::tuple_t<Vs...>, hana::trait<std::is_const>),
                "Collections should contain non-const variable types.");

  constexpr const hana::tuple<Vs...>& variables() const noexcept {
    return *this;
  }

  /// Note: We need this because there is no constexpr lambd a
  struct plus_size {
    template <typename T>
    constexpr std::size_t operator()(std::size_t n, T collection) const
        noexcept {
      return n + collection_traits<T>::size(collection);
    }
  };

  constexpr std::size_t size() const noexcept {
    return hana::fold(variables(), std::size_t{0}, plus_size{});
  }
};

template <typename... Vs>
struct is_collection<static_collection<Vs...>> : bool_constant<true> {};

template <typename... Vs>
constexpr bool operator==(const static_collection<Vs...>& lhs,
                          const static_collection<Vs...>& rhs) noexcept {
  return lhs.variables() == rhs.variables();
}

template <typename... Vs>
constexpr bool operator!=(const static_collection<Vs...>& lhs,
                          const static_collection<Vs...>& rhs) noexcept {
  return !(lhs == rhs);
}

class dynamic_collection;

class dynamic_iterable_variable {
public:
  dynamic_iterable_variable() = default;

  explicit constexpr dynamic_iterable_variable(std::ptrdiff_t value) noexcept
      : m_value{value} {}

  constexpr std::ptrdiff_t value() const noexcept { return m_value; }

private:
  std::ptrdiff_t m_value;
};

#define FUB_DEFINE_CMP_OPERATOR_DYNAMIC_ITERABLE_VARIABLE(op)                  \
  constexpr bool operator op(const dynamic_iterable_variable& lhs,             \
                             const dynamic_iterable_variable& rhs) noexcept {  \
    return lhs.value() op rhs.value();                                         \
  }                                                                            \
  /* */

FUB_DEFINE_CMP_OPERATOR_DYNAMIC_ITERABLE_VARIABLE(<);
FUB_DEFINE_CMP_OPERATOR_DYNAMIC_ITERABLE_VARIABLE(<=);
FUB_DEFINE_CMP_OPERATOR_DYNAMIC_ITERABLE_VARIABLE(>);
FUB_DEFINE_CMP_OPERATOR_DYNAMIC_ITERABLE_VARIABLE(>=);
FUB_DEFINE_CMP_OPERATOR_DYNAMIC_ITERABLE_VARIABLE(==);
FUB_DEFINE_CMP_OPERATOR_DYNAMIC_ITERABLE_VARIABLE(!=);

#undef FUB_DEFINE_CMP_OPERATOR_DYNAMIC_ITERABLE_VARIABLE

class dynamic_collection {
public:
  dynamic_collection() = default;

  explicit constexpr dynamic_collection(std::ptrdiff_t size) noexcept
      : m_size{size} {}

  static constexpr auto variables() noexcept {
    return hana::make_tuple(dynamic_iterable_variable{});
  }

  constexpr std::ptrdiff_t size() const noexcept { return m_size; }

  constexpr dynamic_iterable_variable next(dynamic_iterable_variable var) const
      noexcept {
    return dynamic_iterable_variable{var.value() + 1};
  }

  constexpr dynamic_iterable_variable begin() const noexcept {
    return dynamic_iterable_variable{0};
  }

  constexpr dynamic_iterable_variable end() const noexcept {
    return dynamic_iterable_variable{m_size};
  }

public:
  std::ptrdiff_t m_size;
};

constexpr bool operator==(const dynamic_collection& lhs,
                          const dynamic_collection& rhs) noexcept {
  return lhs.size() == rhs.size();
}

constexpr bool operator!=(const dynamic_collection& lhs,
                          const dynamic_collection& rhs) noexcept {
  return lhs.size() != rhs.size();
}

static_assert(
    std::is_nothrow_copy_constructible<dynamic_iterable_variable>::value,
    "dynamic_iterable_variable is not nothrow copy constructible!");

static_assert(is_nothrow_equality_comparable<dynamic_iterable_variable>::value,
              "dynamic_iterable_variable is not equality comparable!");

static_assert(std::is_nothrow_copy_constructible<dynamic_collection>::value,
              "dynamic_collection is not nothrow copy constructible!");
static_assert(is_nothrow_equality_comparable<dynamic_collection>::value,
              "dynamic_collection is not equality comparable!");

static_assert(
    variable_traits<dynamic_iterable_variable>::is_iterable_in(
        dynamic_collection(1), dynamic_iterable_variable(0)),
    "dynamic_iterable_variable is not itrable in dynamic_collection!");

static_assert(
    collection_traits<dynamic_collection>::size(dynamic_collection{}) == 0,
    "collection size does not match in the trait function.");

static_assert(
    collection_traits<dynamic_collection>::size(dynamic_collection(2)) == 2,
    "collection size does not match in the trait function.");

// }}}

} // namespace fub

#endif
