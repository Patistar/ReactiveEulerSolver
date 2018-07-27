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

#ifndef FUB_GRID_VARIABLE_LIST_HPP
#define FUB_GRID_VARIABLE_LIST_HPP

#include "fub/variable.hpp"

#include <boost/hana/optional.hpp>
#include <boost/hana/range.hpp>
#include <boost/hana/tuple.hpp>
#include <boost/hana/type.hpp>

namespace fub {
inline namespace v1 {
namespace detail {
#ifndef FUB_WITH_CONSTEXPR_LAMBDA
// This is a helper function object to be callable in constexpr contexts.
// It can be replaced in C++17 with a constexpr lambda object.
// see class to_tag_tuple
struct plus_size {
  template <typename V>
  constexpr std::ptrdiff_t operator()(std::ptrdiff_t size, V variable) const
      noexcept {
    return size + variable.size();
  }
};

// This is a helper function object to be callable in constexpr contexts.
// It can be replaced in C++17 with a constexpr lambda object.
// see class to_tag_tuple
struct has_dynamic_size {
  template <typename V>
  constexpr auto operator()(hana::basic_type<V>) const noexcept {
    return hana::integral_constant<bool, V::static_size() == dynamic_extent>{};
  }
};

// This is a helper function object to be callable in constexpr contexts.
// It can be replaced in C++17 with a constexpr lambda object.
// see class to_tag_tuple
struct plus_static_size {
  template <typename V>
  constexpr std::ptrdiff_t operator()(std::ptrdiff_t size,
                                      hana::basic_type<V>) const noexcept {
    return size + V::static_size();
  }
};
#endif // FUB_WITH_CONSTEXPR_LAMBDAS

/// \internal
/// This is a function object which transforms a variable type into a tuple of
/// different valid tag types.
///
/// *Example*:
/// \code
///   struct Density : vector_variable<1> {};
///   struct Momentum : vector_variable<2> {};
///   struct Species : vector_variable<dynamic_extent> {};
///
///   // Returns hana::tuple<tag_t<Density>>
///   to_tag_tuple{}(hana::type_c<Density>);
///
///   // Returns hana::tuple<tag_t<Momentum, 0>, tag_t<Momentum, 1>>
///   to_tag_tuple{}(hana::type_c<Momentum>);
///
///   // Returns hana::tuple<tag_t<Species, dynamic_extent>>
///   to_tag_tuple{}(hana::type_c<Species>);
/// \endcode
class to_tag_tuple {
private:
  template <typename V> struct to_tag {
    template <std::size_t... Is>
    constexpr auto operator()(hana::size_t<Is>...) const noexcept {
      return hana::make_tuple(tag<V, Is>...);
    }
  };

  template <typename V, typename = std::enable_if_t<is_variable<V>{}>>
  static constexpr hana::tuple<tag_t<V, dynamic_extent>>
  construct_tag_tuple(index_constant<dynamic_extent>,
                      hana::basic_type<V>) noexcept {
    return hana::make_tuple(tag<V, dynamic_extent>);
  }

  template <typename V, std::ptrdiff_t N,
            typename = std::enable_if_t<is_variable<V>{}>>
  static constexpr auto construct_tag_tuple(index_constant<N>,
                                            hana::basic_type<V>) noexcept {
    return hana::unpack(hana::make_range(hana::size_c<0>, hana::size_c<N>),
                        to_tag<V>{});
  }

public:
  template <typename V, typename = std::enable_if_t<is_variable<V>{}>>
  constexpr auto operator()(hana::basic_type<V> variable) {
    return construct_tag_tuple(index_c<V::static_size()>, variable);
  }
};

/// This class extents `variant` by defining equality operators with other `tag`
/// types.
template <typename... Tags> struct tag_variant : variant<Tags...> {
  using variant<Tags...>::variant;

  template <typename Tag>
  static constexpr auto is_valid_tag(hana::basic_type<Tag>) noexcept {
    return hana::find(hana::tuple_t<Tags...>, hana::type_c<Tag>) !=
           hana::nothing;
  }

  template <typename Tag>
  static constexpr bool is_valid_tag_v = is_valid_tag(hana::type_c<Tag>);
};

template <typename Tag, typename... Tags,
          typename std::enable_if_t<
              tag_variant<Tags...>::template is_valid_tag_v<Tag>>* = nullptr>
constexpr bool operator==(const Tag& tag,
                          const tag_variant<Tags...>& variant) noexcept {
  if (auto value = get_if<Tag>(&variant)) {
    return *value == tag;
  }
  return false;
}

template <typename... Tags, typename Tag,
          typename std::enable_if_t<
              tag_variant<Tags...>::template is_valid_tag_v<Tag>>* = nullptr>
constexpr bool operator==(const tag_variant<Tags...>& variant,
                          const Tag& tag) noexcept {
  return tag == variant;
}

template <typename... Tags, typename Tag,
          typename std::enable_if_t<
              tag_variant<Tags...>::template is_valid_tag_v<Tag>>* = nullptr>
constexpr bool operator!=(const Tag& tag,
                          const tag_variant<Tags...>& variant) noexcept {
  return !(tag == variant);
}

template <typename... Tags, typename Tag,
          typename std::enable_if_t<
              tag_variant<Tags...>::template is_valid_tag_v<Tag>>* = nullptr>
constexpr bool operator!=(const tag_variant<Tags...>& variant,
                          const Tag& tag) noexcept {
  return !(tag == variant);
}

template <bool, typename... Variables> class variable_list_impl {};

template <typename... Variables>
struct variable_list_impl<true, Variables...> : hana::tuple<Variables...> {
  using Base = hana::tuple<Variables...>;

  template <std::size_t I>
  constexpr std::ptrdiff_t partial_size(hana::size_t<I>) const noexcept {
    return hana::fold_left(hana::take_front(as_tuple(), hana::size_c<I>),
                           std::ptrdiff_t(0), plus_size());
  }

  template <typename Variable, std::ptrdiff_t N>
  constexpr optional<std::ptrdiff_t> index_impl(hana::optional<>,
                                                basic_tag<Variable, N>) const
      noexcept {
    return {};
  }

  template <std::size_t I, typename Variable>
  constexpr optional<std::ptrdiff_t>
  index_impl(hana::optional<hana::size_t<I>> constant_index,
             basic_tag<Variable, dynamic_extent> tag) const noexcept {
    const std::ptrdiff_t index = tag.index();
    if (0 <= index && index < as_tuple()[*constant_index].size()) {
      return partial_size(*constant_index) + index;
    }
    return {};
  }

  template <std::size_t I, typename Variable, std::ptrdiff_t N>
  constexpr optional<std::ptrdiff_t>
  index_impl(hana::optional<hana::size_t<I>>, basic_tag<Variable, N> tag) const
      noexcept {
    return partial_size(hana::size_c<I>) + tag.index();
  }

  static constexpr std::ptrdiff_t static_size_impl(std::true_type) noexcept {
    return dynamic_extent;
  }

  static constexpr std::ptrdiff_t static_size_impl(std::false_type) noexcept {
    return hana::fold_left(hana::tuple_t<Variables...>, std::ptrdiff_t(0),
                           plus_static_size());
  }

  static constexpr auto as_tags() noexcept {
    return hana::flatten(
        hana::transform(hana::tuple_t<Variables...>, to_tag_tuple()));
  }

  using variant_type =
      list_cast_t<tag_variant,
                  remove_cvref_t<decltype(variable_list_impl::as_tags())>>;

  template <typename V>
  std::tuple<basic_tag<V, dynamic_extent>, std::ptrdiff_t, std::ptrdiff_t>
  increase_index(basic_tag<V, dynamic_extent> tag, std::ptrdiff_t dest_index,
                 std::ptrdiff_t n) const noexcept {
    const std::ptrdiff_t dynamic_size = std::min(n, variable(tag).size());
    return std::make_tuple(basic_tag<V, dynamic_extent>{dynamic_size},
                           dest_index +
                               std::min(n - dynamic_size, std::ptrdiff_t(1)),
                           n - dynamic_size);
  }

  template <typename V, std::ptrdiff_t N>
  std::tuple<basic_tag<V, N>, std::ptrdiff_t, std::ptrdiff_t>
  increase_index(basic_tag<V, N> tag, std::ptrdiff_t dest_index,
                 std::ptrdiff_t n) const noexcept {
    const std::ptrdiff_t step_size = std::min(n, std::ptrdiff_t(1));
    return std::make_tuple(tag, dest_index + step_size, n - step_size);
  }

  template <typename Tag>
  static constexpr auto is_valid_tag(hana::basic_type<Tag>) noexcept {
    return hana::find(as_tags(), hana::type_c<Tag>) != hana::nothing;
  }

  using Base::Base;

  template <typename Tag>
  static constexpr bool is_valid_tag_v = is_valid_tag(hana::type_c<Tag>);

  static constexpr hana::tuple<hana::type<Variables>...>
  as_type_tuple() noexcept {
    return {};
  }

  constexpr hana::tuple<Variables...> as_tuple() const noexcept {
    return *this;
  }

  static constexpr std::ptrdiff_t static_size() noexcept {
    auto any_dynamic_extent =
        hana::any_of(hana::tuple_t<Variables...>, has_dynamic_size());
    return static_size_impl(any_dynamic_extent);
  }

  constexpr std::ptrdiff_t size() const noexcept {
    return partial_size(hana::size(as_tuple()));
  }

  template <typename Variable, std::ptrdiff_t N>
  constexpr optional<std::ptrdiff_t> index(basic_tag<Variable, N> tag) const
      noexcept {
    return index_impl(hana::index_if(to_types(as_tuple()),
                                     hana::equal.to(hana::type_c<Variable>)),
                      tag);
  }

  template <typename Variable, std::ptrdiff_t N>
  constexpr Variable variable(basic_tag<Variable, N>) const noexcept {
    auto index = hana::index_if(hana::tuple_t<Variables...>,
                                hana::equal.to(hana::type_c<Variable>));
    return as_tuple()[*index];
  }

  constexpr optional<variant_type> tag(std::ptrdiff_t n) const noexcept {
    // GCC 6.3 needs to explicitly call this with variable_list_impl::as_tags
    constexpr auto tags = variable_list_impl::as_tags();
    constexpr auto size = hana::size(tags);
    variant_type tag_array[size];
    std::ptrdiff_t counter = n;
    std::ptrdiff_t dest_index = 0;
    hana::for_each(hana::make_range(hana::size_c<0>, size),
                   [&, this](auto index) {
                     std::tie(tag_array[index], dest_index, counter) =
                         this->increase_index(tags[index], dest_index, counter);
                   });
    if (0 <= dest_index && dest_index < size()) {
      return tag_array[dest_index];
    }
    return {};
  }
};
} // namespace detail

/// A variable_list is a searchable collection of variables.
///
/// Whenever the size of all Variable are compile-time known this class provides
/// a stateless and compile-time enumeration of these variables. This class
/// proviedes methods to inspect and visit all the variables which it contains.
///
/// \example variable_list.cpp
template <typename... Variables>
class variable_list
    : private detail::variable_list_impl<
          conjunction<is_variable<Variables>...>::value, Variables...> {
  using Base = typename detail::variable_list_impl<
      conjunction<is_variable<Variables>...>::value, Variables...>;

public:
  using variant_type = typename Base::variant_type;

  using Base::Base;

  /// \name Static Observers

  /// This is true if Tag is a valid tag type to be used with this variable
  /// list.
  template <typename Tag>
  static constexpr bool is_valid_tag_v = is_valid_tag(hana::type_c<Tag>);

  /// \return Returns a tuple of type objects.
  ///
  /// \throws Nothing.
  static constexpr hana::tuple<hana::type<Variables>...>
  as_type_tuple() noexcept {
    return {};
  }

  /// Returns the compile-time known number of all indexable variable elements.
  ///
  /// \post static_size() > 0 || static_size() == dynamic_extent.
  ///
  /// \throws Nothing.
  static constexpr std::ptrdiff_t static_size() noexcept {
    return Base::static_size();
  }

  /// \name Observers

  /// \return Returns a tuple of contained variables.
  ///
  /// This function will return all variables with its correctly set sizes, even
  /// if they are set in run-time only.
  ///
  /// \throws Nothing.
  constexpr hana::tuple<Variables...> as_tuple() const noexcept {
    return Base::as_tuple();
  }

  /// \return Returns the total accumulated size of all variables.
  ///
  /// This is also the codomain size of index(tag<Varbiable>), i.e. it is equal
  /// to the number of valid tags which can be used to reference to the vairable
  /// elements.
  ///
  /// \throws Nothing.
  constexpr std::ptrdiff_t size() const noexcept { return Base::size(); }

  template <typename Variable, std::ptrdiff_t N>
  constexpr optional<std::ptrdiff_t> index(basic_tag<Variable, N> tag) const
      noexcept {
    return Base::index(tag);
  }

  template <typename Variable, std::ptrdiff_t N>
  constexpr Variable variable(basic_tag<Variable, N> tag) const noexcept {
    return Base::variable(tag);
  }

  constexpr optional<variant_type> tag(std::ptrdiff_t n) const noexcept {
    return Base::tag(n);
  }
};

template <typename... Variables>
constexpr bool operator==(const variable_list<Variables...>& lhs,
                          const variable_list<Variables...>& rhs) noexcept {
  return lhs.as_tuple() == rhs.as_tuple();
}

namespace detail {
template <typename L> using as_tuple_t = decltype(std::declval<L>().as_tuple());
}

/// \ingroup type-traits
/// Is true_type if `L` fulfills the \b VariableList concept.
template <typename L>
struct is_variable_list
    : conjunction<is_variable<L>, is_detected<detail::as_tuple_t, L>> {};

template <typename... Variables>
constexpr bool operator!=(const variable_list<Variables...>& lhs,
                          const variable_list<Variables...>& rhs) noexcept {
  return lhs.as_tuple() != rhs.as_tuple();
}

namespace detail {
template <typename F, typename Variable> struct invoke_function_with_tag {
  F* function;
  template <typename I> constexpr void operator()(I index) const noexcept {
    fub::invoke(*function, tag<Variable, index>);
  }
};

template <std::ptrdiff_t N, typename V, typename F,
          typename = std::enable_if_t<is_variable<V>::value>>
constexpr F for_each_helper(index_constant<N>, V, F function) {
  hana::for_each(hana::make_range(hana::size_c<0>, hana::size_c<N>),
                 invoke_function_with_tag<F, V>{&function});
  return function;
}

template <typename V, typename F,
          typename = std::enable_if_t<is_variable<V>::value>>
constexpr F for_each_helper(index_constant<dynamic_extent>, V variable,
                            F function) {
  const std::ptrdiff_t size = variable.size();
  for (std::ptrdiff_t i = 0; i < size; ++i) {
    fub::invoke(function, tag_t<V>(i));
  }
  return function;
}
} // namespace detail

/// This overload iterates through every distinct variable tag and invokes
/// `function(tag)`.
template <typename V, typename F,
          typename = std::enable_if_t<!is_variable_list<V>::value>,
          typename = std::enable_if_t<is_variable<V>::value>>
constexpr F for_each(V variable, F function) {
  return detail::for_each_helper(index_c<V::static_size()>, variable, function);
}

/// This overload iterates through every distinct variable tag for all contained
/// variables in the variable_list.
template <typename List, typename F,
          typename = std::enable_if_t<is_variable_list<List>{}>>
F for_each(const List& list, F function) {
  hana::for_each(list.as_tuple(), [&](auto variable) {
    for_each(variable, std::ref(function));
  });
  return function;
}

} // namespace v1
} // namespace fub

#endif