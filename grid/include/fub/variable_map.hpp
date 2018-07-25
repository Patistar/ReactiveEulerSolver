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
#include "fub/variant.hpp"

#include <boost/hana.hpp>

namespace fub {

namespace hana = boost::hana;

template <std::ptrdiff_t Rank> class vector_variable : private extents<Rank> {
private:
  using Base = extents<Rank>;

public:
  using Base::Base;

  static constexpr std::ptrdiff_t static_size() noexcept {
    return Base::static_extent(0);
  }

  constexpr std::ptrdiff_t size() const noexcept { return Base::extent(0); }
};
using scalar_variable = vector_variable<1>;

template <std::ptrdiff_t RankL, std::ptrdiff_t RankR>
constexpr bool operator==(const vector_variable<RankL>& lhs,
                          const vector_variable<RankR>& rhs) noexcept {
  return lhs.size() == rhs.size();
}

template <std::ptrdiff_t RankL, std::ptrdiff_t RankR>
constexpr bool operator!=(const vector_variable<RankL>& lhs,
                          const vector_variable<RankR>& rhs) noexcept {
  return lhs.size() != rhs.size();
}

namespace detail {

template <typename T> using static_size_t = decltype(T::static_size());
template <typename T> using size_t = decltype(std::declval<T>().size());
template <typename T>
struct is_variable : conjunction<is_detected<static_size_t, T>,
                                 is_detected<size_t, const T&>> {};

template <bool IsVariable, typename V, std::ptrdiff_t N>
class tag_check_variable {
  static_assert(IsVariable, "V does not fulfill the variable concept.");
};

template <bool OutOfRange, typename V, std::ptrdiff_t N>
class tag_check_out_of_range {
  static_assert(OutOfRange,
                "N is out of range, i.e. dimension of V is lower than N.");
};

template <std::ptrdiff_t N, std::ptrdiff_t Size>
struct is_both_dynamic_or_less_than
    : std::integral_constant<
          bool, ((N == dynamic_extent && Size == dynamic_extent) || N < Size)> {
};

template <typename V, std::ptrdiff_t N>
class tag_check_variable<true, V, N>
    : public tag_check_out_of_range<
          is_both_dynamic_or_less_than<N, V::static_size()>::value, V, N> {
  using tag_check_out_of_range<
      is_both_dynamic_or_less_than<N, V::static_size()>::value, V,
      N>::tag_check_out_of_range;
};

template <typename V, std::ptrdiff_t N>
class tag_check_out_of_range<true, V, N> : private extents<N> {
public:
  using extents<N>::extents;

  static constexpr std::ptrdiff_t static_index() noexcept {
    return extents<N>::static_extent(0);
  }
  constexpr std::ptrdiff_t index() const noexcept {
    return extents<N>::extent(0);
  }
};

template <typename V, std::ptrdiff_t N>
class basic_tag : public tag_check_variable<is_variable<V>::value, V, N> {
  using tag_check_variable<is_variable<V>::value, V, N>::tag_check_variable;
};

template <typename V, std::ptrdiff_t N = 0>
using tag_t = std::conditional_t<V::static_size() == dynamic_extent,
                                 basic_tag<V, dynamic_extent>, basic_tag<V, N>>;

template <typename T, std::ptrdiff_t N, typename S, std::ptrdiff_t M>
constexpr bool operator==(const basic_tag<T, N>& lhs,
                          const basic_tag<S, M>& rhs) noexcept {
  return lhs.index() == rhs.index() && std::is_same<T, S>{};
}

template <typename T, std::ptrdiff_t N, typename S, std::ptrdiff_t M>
constexpr bool operator!=(const basic_tag<T, N>& lhs,
                          const basic_tag<S, M>& rhs) noexcept {
  return !(lhs == rhs);
}

template <typename V, std::ptrdiff_t N = 0> static constexpr tag_t<V, N> tag{};

struct plus_size {
  template <typename V>
  constexpr std::ptrdiff_t operator()(std::ptrdiff_t size, V variable) const
      noexcept {
    return size + variable.size();
  }
};

struct is_dynamic_extent {
  template <typename V>
  constexpr auto operator()(hana::basic_type<V>) const noexcept {
    return hana::integral_constant<bool, V::static_size() == dynamic_extent>{};
  }
};

struct plus_static_size {
  template <typename V>
  constexpr std::ptrdiff_t operator()(std::ptrdiff_t size,
                                      hana::basic_type<V>) const noexcept {
    return size + V::static_size();
  }
};

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

template <typename... Tags> struct tag_variant : public variant<Tags...> {
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
    if (auto value = std::get_if<Tag>(&variant)) {
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

  using variant_type = list_cast_t<tag_variant, decltype(as_tags())>;

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
  constexpr auto variable(basic_tag<Variable, N>) const noexcept {
    auto index = hana::index_if(hana::tuple_t<Variables...>,
                                hana::equal.to(hana::type_c<Variable>));
    return as_tuple()[*index];
  }

  constexpr optional<variant_type> variable(std::ptrdiff_t n) const noexcept {
    constexpr auto tags = as_tags();
    constexpr auto size = hana::size(tags);
    variant_type tag_array[size];
    std::ptrdiff_t counter = n;
    std::ptrdiff_t dest_index = 0;
    hana::for_each(hana::make_range(hana::size_c<0>, size), [&](auto index) {
      std::tie(tag_array[index], dest_index, counter) =
          increase_index(tags[index], dest_index, counter);
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

template <typename Span, typename Mapping>
class basic_variable_mdspan_storage : Mapping {
public:
  basic_variable_mdspan_storage() = default;
  basic_variable_mdspan_storage(Span span, Mapping mapping) noexcept
      : Mapping(mapping), m_span(span) {}

  constexpr Span span() const noexcept { return m_span; }

  constexpr Mapping get_mapping() const noexcept { return *this; }

private:
  Span m_span;
};

template <typename Mapping, typename VariableList> struct compute_span_extents {
  static constexpr std::ptrdiff_t
  extent(hana::basic_type<Mapping>, hana::basic_type<VariableList>) noexcept {
    if (Mapping::static_required_span_size() == dynamic_extent ||
        VariableList::static_size() == dynamic_extent) {
      return dynamic_extent;
    }
    return Mapping::static_required_span_size() * VariableList::static_size();
  }

  using type =
      extents<extent(hana::type_c<Mapping>, hana::type_c<VariableList>)>;
};

template <typename Mapping, typename VariableList>
using compute_span_extents_t =
    typename compute_span_extents<Mapping, VariableList>::type;

template <typename VariableList, typename T, typename Extents,
          typename LayoutPolicy, typename Accessor>
class basic_variable_mdspan : VariableList {
public:
  using variable_list = VariableList;
  using extents = Extents;
  using accessor = Accessor;
  using layout_policy = LayoutPolicy;
  using mapping = typename LayoutPolicy::template mapping<Extents>;
  using element_type = typename Accessor::element_type;
  using pointer = typename Accessor::pointer;
  using mdspan = basic_mdspan<element_type, extents, layout_policy, accessor>;
  using span_extents = compute_span_extents_t<mapping, variable_list>;
  using span_type =
      basic_span<element_type, span_extents::static_extent(0), accessor>;

  template <typename Tag>
  static constexpr bool is_valid_tag_v =
      variable_list::template is_valid_tag_v<Tag>;

  basic_variable_mdspan() = default;

  constexpr basic_variable_mdspan(span_type data, variable_list list,
                                  extents e = extents()) noexcept
      : VariableList(list), m_storage(data, mapping(e)) {
    assert(span().size() >= required_span_size());
  }

  template <typename Variable, std::ptrdiff_t N,
            typename std::enable_if_t<is_valid_tag_v<Variable>>* = nullptr>
  constexpr mdspan operator[](basic_tag<Variable, N> tag) const noexcept {
    const variable_list variables = get_variable_list();
    optional<std::ptrdiff_t> index = variables.index(tag);
    assert(index);
    mapping map = get_mapping();
    const std::ptrdiff_t size = map.required_span_size();
    return mdspan(subspan(span(), *index * size, size), map);
  }

  constexpr variable_list get_variable_list() const noexcept { return *this; }

  constexpr span_type span() const noexcept { return m_storage.span(); }

  constexpr mapping get_mapping() const noexcept {
    return m_storage.get_mapping();
  }

  constexpr std::ptrdiff_t required_span_size() const noexcept {
    return get_variable_list().size() * get_mapping().required_span_size();
  }

  friend constexpr bool operator==(const basic_variable_mdspan& lhs,
                                   const basic_variable_mdspan& rhs) noexcept {
    return lhs.span() == rhs.span() &&
           lhs.get_variable_list() == rhs.get_variable_list() &&
           lhs.get_mapping() == rhs.get_mapping();
  }

  friend constexpr bool operator!=(const basic_variable_mdspan& lhs,
                                   const basic_variable_mdspan& rhs) noexcept {
    return !(lhs == rhs);
  }

private:
  basic_variable_mdspan_storage<span_type, mapping> m_storage;
};

template <typename VariableList, typename T, typename Extents>
using variable_mdspan = basic_variable_mdspan<VariableList, T, Extents,
                                              layout_left, accessor_native<T>>;

} // namespace detail

using detail::basic_tag;
using detail::for_each;
using detail::is_variable;
using detail::tag;
using detail::tag_t;
using detail::variable_list;
using detail::basic_variable_mdspan;

} // namespace fub