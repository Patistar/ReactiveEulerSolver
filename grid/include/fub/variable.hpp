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

#include "fub/extents.hpp"
#include "fub/functional.hpp"
#include "fub/hana.hpp"
#include "fub/mdspan.hpp"
#include "fub/optional.hpp"
#include "fub/type_traits.hpp"
#include "fub/variant.hpp"

#include <boost/hana.hpp>

namespace fub {
namespace v1 {
namespace hana = boost::hana;

template <typename T> using static_size_fn = decltype(T::static_size());
template <typename T> using size_memfn = decltype(std::declval<T>().size());

/// This is a trait type which checks if the `T` fulfills the \b Variable
/// requirements.
///
/// \see vector_variable
template <typename T>
struct is_variable : conjunction<is_regular<T>, is_detected<static_size_fn, T>,
                                 is_detected<size_memfn, const T&>> {};

/// This is a helper class which make it easier to create types which satisfy
/// the \b Variable concept.
///
/// The intent is to derive from this class and to automatically satisfy all
/// requirements to be used correctly within this framework.
///
/// *Example:*
///
/// \code
/// // Density is a scalar variable
/// struct Density : vector_variable<1> {};
///
/// // Momentum is vector variable of rank three, i.e. is has three elements.
/// struct Momentum : vector_variable<3> {};
///
/// // Species is a vector variable of dynamic rank. i.e. its size is determined
/// // at runtime.
/// struct Species : vector_variable<dynamic_extent> {};
/// \endcode
///
/// \tparam Rank The size of the variable.
///
/// The template parameter Rank can be either a positive integral or
/// `dynamic_extent` which is for variables which have their size determined at
/// run-time.
template <std::ptrdiff_t Rank> class vector_variable : private extents<Rank> {
public:
  using extents<Rank>::extents;

  /// Returns the compile-time Rank of the variable.
  ///
  /// This is usefull in generic contexts if you want to check if a given
  /// variable is statically sized at compile-time or dynamically sized at
  /// run-time.
  ///
  /// \throws Nothing.
  static constexpr std::ptrdiff_t static_size() noexcept {
    return Base::static_extent(0);
  }

  /// Returns the run-time Rank of the variable.
  ///
  /// \throws Nothing.
  constexpr std::ptrdiff_t size() const noexcept { return Base::extent(0); }
};

/// This is a conventient typedef for one-dimensional vector variables.
using scalar_variable = vector_variable<1>;

/// Returns true if two vector variables have the same size.
///
/// \note vector_variables are immutable, thus it makes sense to compare their
/// run-time size and not their compile-time size.
template <std::ptrdiff_t RankL, std::ptrdiff_t RankR>
constexpr bool operator==(const vector_variable<RankL>& lhs,
                          const vector_variable<RankR>& rhs) noexcept {
  return lhs.size() == rhs.size();
}

/// Returns true if two vector variables have different sizes.
template <std::ptrdiff_t RankL, std::ptrdiff_t RankR>
constexpr bool operator!=(const vector_variable<RankL>& lhs,
                          const vector_variable<RankR>& rhs) noexcept {
  return lhs.size() != rhs.size();
}

// This trait checks if `N` is a valid size
template <std::ptrdiff_t N, std::ptrdiff_t Size>
struct is_both_dynamic_or_less_than
    : std::integral_constant<
          bool, ((N == dynamic_extent && Size == dynamic_extent) || N < Size)> {
};

/// This type is used to index into elements of a `V`.
///
/// \tparam V needs to fulfill the Variable concept.
/// \tparam N the element index of the element in `V`.
///
/// N can be either a positive integral or `dynamic_extent`.
template <typename V, std::ptrdiff_t N> class basic_tag : extents<N> {
public:
  using extents<N>::extents;

  static_assert(is_variable<V>::value,
                "V does not fulfill the Variable concept.");

  static_assert(is_both_dynamic_or_less_than<N, V::static_size()>::value,
                "N is out of range.");

  static constexpr std::ptrdiff_t static_index() noexcept {
    return extents<N>::static_extent(0);
  }
  constexpr std::ptrdiff_t index() const noexcept {
    return extents<N>::extent(0);
  }
};

/// This is a short-hand automatically create dynamic-sized tags from
/// dynamic-sized variables.
template <typename V, std::ptrdiff_t N = 0>
using tag_t = std::conditional_t<V::static_size() == dynamic_extent,
                                 basic_tag<V, dynamic_extent>, basic_tag<V, N>>;

/// Tags are equal if they have the same index and variable.
template <typename T, std::ptrdiff_t N, typename S, std::ptrdiff_t M>
constexpr bool operator==(const basic_tag<T, N>& lhs,
                          const basic_tag<S, M>& rhs) noexcept {
  return lhs.index() == rhs.index() && std::is_same<T, S>{};
}

/// This is equivalent to `!(lhs == rhs)`.
template <typename T, std::ptrdiff_t N, typename S, std::ptrdiff_t M>
constexpr bool operator!=(const basic_tag<T, N>& lhs,
                          const basic_tag<S, M>& rhs) noexcept {
  return !(lhs == rhs);
}

template <typename V, std::ptrdiff_t N = 0> static constexpr tag_t<V, N> tag{};

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
struct is_dynamically_sized {
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

  template <typename Tag, typename = std::enable_if_t<is_valid_tag_v<Tag>>>
  friend constexpr bool operator==(const Tag& tag,
                                   const tag_variant& variant) noexcept {
    if (auto value = get_if<Tag>(&variant)) {
      return *value == tag;
    }
    return false;
  }

  template <typename Tag, typename = std::enable_if_t<is_valid_tag_v<Tag>>>
  friend constexpr bool operator==(const tag_variant& variant,
                                   const Tag& tag) noexcept {
    return tag == variant;
  }

  template <typename Tag, typename = std::enable_if_t<is_valid_tag_v<Tag>>>
  friend constexpr bool operator!=(const Tag& tag,
                                   const tag_variant& variant) noexcept {
    return !(tag == variant);
  }

  template <typename Tag, typename = std::enable_if_t<is_valid_tag_v<Tag>>>
  friend constexpr bool operator!=(const tag_variant& variant,
                                   const Tag& tag) noexcept {
    return !(tag == variant);
  }
};

template <bool, typename... Variables> class variable_list_impl {};

template <typename... Variables>
class variable_list_impl<true, Variables...>
    : private hana::tuple<Variables...> {
private:
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

public:
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
        hana::any_of(hana::tuple_t<Variables...>, is_dynamic_extent());
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

  constexpr optional<variant_type> variable(std::ptrdiff_t n) const noexcept {
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

template <typename... Variables>
class variable_list
    : public variable_list_impl<conjunction<is_variable<Variables>...>::value,
                                Variables...> {
  using variable_list_impl<conjunction<is_variable<Variables>...>::value,
                           Variables...>::variable_list_impl;
};

template <typename... Variables>
constexpr bool operator==(const variable_list<Variables...>& lhs,
                          const variable_list<Variables...>& rhs) noexcept {
  return lhs.as_tuple() == rhs.as_tuple();
}

template <typename... Variables>
constexpr bool operator!=(const variable_list<Variables...>& lhs,
                          const variable_list<Variables...>& rhs) noexcept {
  return lhs.as_tuple() != rhs.as_tuple();
}

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

template <typename V, typename F,
          typename = std::enable_if_t<is_variable<V>::value>>
constexpr F for_each(V variable, F function) {
  return for_each_helper(index_c<V::static_size()>, variable, function);
}

template <typename F, typename... Variables>
F for_each(const variable_list<Variables...>& list, F function) {
  hana::for_each(list.as_tuple(), [&](auto variable) {
    for_each(variable, std::ref(function));
  });
  return function;
}

} // namespace v1

using v1::basic_tag;
using v1::basic_variable_mdspan;
using v1::for_each;
using v1::is_variable;
using v1::tag;
using v1::tag_t;
using v1::variable_list;
using v1::variable_mdspan;

} // namespace fub