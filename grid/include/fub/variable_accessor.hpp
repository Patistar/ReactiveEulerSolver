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

#include "fub/accessor_native.hpp"
#include "fub/hana.hpp"
#include "fub/span.hpp"
#include "fub/variables.hpp"

#include <boost/hana/equal.hpp>
#include <boost/hana/ext/std/integral_constant.hpp>
#include <boost/hana/find.hpp>
#include <boost/hana/flatten.hpp>
#include <boost/hana/if.hpp>
#include <boost/hana/index_if.hpp>
#include <boost/hana/not_equal.hpp>
#include <boost/hana/optional.hpp>
#include <boost/hana/or.hpp>
#include <boost/hana/transform.hpp>
#include <boost/hana/tuple.hpp>
#include <boost/hana/type.hpp>

#include <array>

namespace fub {
namespace detail {
template <typename Variables, typename ProtoAccessor, typename DefaultType>
struct variable_accessor_implementation;
}

template <typename Variables, typename ProtoAccessor = accessor_native<>,
          typename DefaultType = double>
class variable_accessor
    : private detail::variable_accessor_implementation<Variables, ProtoAccessor,
                                                       DefaultType> {
private:
  using base =
      detail::variable_accessor_implementation<Variables, ProtoAccessor,
                                               DefaultType>;

public:
  variable_accessor() = default;
  variable_accessor(const variable_accessor&) = default;
  variable_accessor& operator=(const variable_accessor&) = default;
  variable_accessor(variable_accessor&&) = default;
  variable_accessor& operator=(variable_accessor&&) = default;

  template <typename OtherVariablesList>
  using rebind =
      variable_accessor<OtherVariablesList, ProtoAccessor, DefaultType>;

  template <typename OtherVars>
  variable_accessor(
      const variable_accessor<OtherVars, ProtoAccessor, DefaultType>& other) {}

  // Variable dependent typedefs
  template <typename V> using accessor = typename base::template accessor<V>;

  template <typename V>
  using value_type = typename base::template value_type<V>;

  template <typename V>
  using element_type = typename base::template element_type<V>;

  template <typename V> using pointer = typename base::template pointer<V>;

  template <typename V>
  using const_pointer = typename base::template const_pointer<V>;

  template <typename V> using reference = typename base::template reference<V>;

  template <typename V>
  using const_reference = typename base::template const_reference<V>;

  template <typename V> using span = typename base::template span<V>;

  template <typename V>
  using const_span = typename base::template const_span<V>;

  // Storage related typedefs

  using value_storage = typename base::value_storage;
  using pointer_storage = typename base::pointer_storage;
  using const_pointer_storage = typename base::const_pointer_storage;

  /// Returns a tuple containing all types which can be used to access data via
  /// this accessor.
  using base::get_accessible_variables;

  /// Transforms a const value storage tuple into pointer to const storage tuple
  /// where each tuple entry points to the value tuple entry.
  using base::to_pointer_storage;

  using base::make_span;

  /// Returns a reference to a variable value with specified storage, total
  /// extents and a given offset.
  using base::access;

  template <typename V,
            typename = std::enable_if_t<
                base::find_collection(hana::type_c<V>) != hana::nothing>>
  static constexpr auto index(hana::basic_type<V>) noexcept {
    return *base::find_collection_index(hana::type_c<V>);
  }
};

///////////////////////////////////////////////////////////////////////////////
//                                                               Implementation

namespace detail {

template <typename T> struct collection_comparator {
  T m_value;

  template <typename S>
  constexpr bool equal_to(std::true_type,
                          const collection_comparator<S>& rhs) const noexcept {
    return m_value == rhs.m_value;
  }

  template <typename S>
  constexpr bool equal_to(std::false_type,
                          const collection_comparator<S>& rhs) const noexcept {
    return false;
  }

  template <typename S>
  constexpr bool operator==(const collection_comparator<S>& rhs) const noexcept {
    return equal_to(hana::type_c<T> == hana::type_c<S>, rhs);
  }

  template <typename S>
  constexpr bool operator!=(const collection_comparator<S>& rhs) const noexcept {
    return !(*this == rhs);
  }
};

struct make_comparator_t {
  template <typename T>
  constexpr collection_comparator<T> operator()(T value) const noexcept {
    return {value};
  }
};
static constexpr make_comparator_t make_comparator{};

/// Returns all available unions as a hana::tuple object.
///
/// Variables are assumed to be light-weight handles to some data. Thus we
/// require some exception guarantees:
///
/// Require: is_nothrow_copy_constructible<Variable>
/// Require: is_nothrow_equality_comparable<const Variable&>
/// @{
template <typename Variable, bool IsCollection>
struct variable_collections_impl : hana::tuple<Variable> {
  using Base = hana::tuple<Variable>;
  using Base::Base;

  static_assert(std::is_nothrow_copy_constructible<Base>{},
                "hana::tuple<Variable> must be nothrow copy constructible.");

  static_assert(is_equality_comparable<Base>{},
                "hana::tuple<Variable> must be nothrow equality comparable.");

  constexpr hana::tuple<Variable> get_collections() const noexcept {
    return *this;
  }
};

template <typename Pred> struct plus_one_until_t {
  Pred pred;
  using result_type = std::pair<std::ptrdiff_t, bool>;

  template <typename T>
  constexpr result_type operator()(result_type result, T x) const noexcept {
    std::ptrdiff_t n = result.first;
    bool found = result.second || pred(x);
    if (!found) {
      n += 1;
    }
    return result_type(n, found);
  }
};

template <typename Pred>
constexpr plus_one_until_t<Pred> plus_one_until(Pred pred) noexcept {
  return {pred};
}

template <typename Xs, typename Pred>
constexpr std::ptrdiff_t index_if(Xs xs, Pred pred) {
  auto x = hana::front(xs);
  auto rest = hana::drop_front(xs, int_c<1>);
  return hana::fold_left(rest, std::make_pair(std::ptrdiff_t{0}, pred(x)),
                         plus_one_until(pred))
      .first;
}

template <template <typename...> class List, typename... Vs>
struct variable_collections_impl<List<Vs...>, false> : hana::tuple<Vs...> {
  using Base = hana::tuple<Vs...>;
  using Base::Base;

  static_assert(std::is_nothrow_copy_constructible<Base>{},
                "hana::tuple<Vs...> must be nothrow copy constructible.");

  static_assert(is_equality_comparable<Base>{},
                "hana::tuple<Vs...> must be nothrow equality comparable.");

  constexpr hana::tuple<Vs...> get_collections() const noexcept {
    return *this;
  }
};

template <typename Variable>
struct variable_collections : variable_collections_impl<Variable, false> {
  using Base = variable_collections_impl<Variable, false>;
  using Base::Base;
};

template <template <typename...> class List, typename... Vs>
struct variable_collections<List<Vs...>>
    : variable_collections_impl<List<Vs...>,
                                is_collection<List<Vs...>>::value> {
  using Base =
      variable_collections_impl<List<Vs...>, is_collection<List<Vs...>>::value>;
  using Base::Base;
};
/// @}

////////////////////////////////////////////////////////////////////////////////
//                                                                 variable_find

/// Given a Variable which is going to be accessed, this class implements how to
/// find the union where the value is stored at. It also computes the index of
/// the union in the set of all present unions.
///
/// Note: This class template works for both cases: For tuple-like types
/// Variables and simple unions.
template <typename Collections>
struct variable_find : variable_collections<Collections> {
  using base = variable_collections<Collections>;
  using base::base;

  // Define Type based Observers

  /// Returns a tuple of types containing all collection types accessible by
  /// this accessor.
  ///
  /// Throws: Nothing.
  static constexpr auto get_collection_types() noexcept {
    return decltype(to_types(
        std::declval<variable_collections<Collections>>().get_collections())){};
  }

  /// Returns an index to the specified collection `C` if it exists in
  /// `Collections`. If it does not exist the function returns `hana::nothing`.
  ///
  /// Throws: Nothing.
  template <typename C>
  static constexpr auto index_collection_type(hana::basic_type<C> c) noexcept {
    return hana::index_if(get_collection_types(), hana::equal.to(c));
  }

  /// This is a unary predicate which checks if a variable is accessible in a
  /// collection.
  template <typename Variable> struct is_accessible_type {
    template <typename Collection>
    constexpr auto operator()(hana::basic_type<Collection> collection) {
      hana::type<Variable> variable;
      return hana::or_(
          collection == variable,
          hana::find(
              to_types(collection_traits<Collection>::variables(Collection{})),
              variable) != hana::nothing);
    }
  };

  /// Returns an optional index if the specified variable type is accessible for
  /// the collection type referred by the returned index.
  ///
  /// Throws: Nothing.
  template <typename Variable>
  static constexpr auto
  index_accessible_collection_type(hana::basic_type<Variable> variable =
                                       hana::basic_type<Variable>()) noexcept {
    return hana::index_if(variable_find::get_collection_types(),
                          variable_find::is_accessible_type<Variable>{});
  }

  template <typename Collection>
  constexpr auto find_collection(Collection collection) const noexcept {
    return hana::find(base::get_collections(), collection);
  }

  /// This is a unary predicate type which checks if a variable is accessible in
  /// a collection.
  template <typename Variable> struct is_accessible {
    collection_comparator<Variable> variable;

    template <typename Collection>
    constexpr auto operator()(Collection collection) const noexcept {
      return hana::or_(
          collection_comparator<Collection>{collection} == variable,
          hana::index_if(
              hana::transform(
                  collection_traits<Collection>::variables(collection),
                  make_comparator),
              hana::equal.to(variable)) != hana::nothing);
    }
  };

  /// Returns an optional index `I` for which `V` is accessible in
  /// `get_collections()[*I]`.
  ///
  /// Throws: Nothing.
  template <typename Variable>
  constexpr auto index_accessible_collection(Variable variable) const noexcept {
    return index_if(base::get_collections(), is_accessible<Variable>{variable});
  }

  /// Returns an optional union varibale `u` for which `variable` is accessible.
  ///
  /// Throws: Nothing.
  template <typename Variable>
  constexpr auto find_accessible_collection(Variable variable) const noexcept {
    return hana::find_if(base::get_collections(),
                         is_accessible<Variable>{variable});
  }

  template <typename V,
            typename std::enable_if_t<index_accessible_collection_type<V>() !=
                                      hana::nothing>* = nullptr>
  constexpr auto local_index(hana::basic_type<V> var) const noexcept {
    constexpr auto remove_const = hana::metafunction<std::remove_const>;
    auto index = index_accessible_collection(var);
    using Collection = decltype(base::get_collections()[*index]);
    return hana::index_if(collection_traits<Collection>::variables(),
                          hana::equal.to(remove_const(var)));
  }

  /// Returns a tuple of all accessible variables.
  ///
  /// Throws: Nothing.
  constexpr auto get_variables() const noexcept {
    return hana::flatten(
        hana::transform(base::get_collections(), [](auto collection) {
          using C = remove_cvref_t<decltype(collection)>;
          return collection_traits<C>::variables(collection);
        }));
  }

  /// Returns a tuple of all accessible variable types.
  ///
  /// Throws: Nothing.
  constexpr auto get_variable_types() const noexcept {
    return hana::flatten(
        hana::transform(get_collection_types(), [](auto collection) {
          using C = remove_cvref_t<typename decltype(collection)::type>;
          return decltype(
              to_types(collection_traits<C>::variables(std::declval<C>()))){};
        }));
  }
};

///////////////////////////////////////////////////////////////////////////////
//                                                   variable_accessor_typedefs
//
// Base case: VariableUnion is a single union of variables

template <typename VariableUnion, typename ProtoAccessor, typename DefaultType>
struct variable_accessor_typedefs : variable_find<VariableUnion> {
  using base = variable_find<VariableUnion>;
  using base::contains_collection;
  using base::find_collection;
  using base::find_collection_index;

  template <typename V> using value_type_t = typename V::value_type;

  using variable_value_type =
      detected_or_t<DefaultType, value_type_t, VariableUnion>;

  template <typename V>
  using variable_element_type = std::conditional_t<
      std::is_const<V>::value || std::is_const<VariableUnion>::value,
      std::add_const_t<variable_value_type>, variable_value_type>;

  template <
      typename V,
      typename std::enable_if_t<disjunction<
          bool_constant<find_collection(hana::type_c<V>) != hana::nothing>,
          std::is_same<remove_cvref_t<V>,
                       remove_cvref_t<VariableUnion>>>::value>* = nullptr>
  static constexpr auto get_accessor(hana::basic_type<V>) {
    using Result =
        typename ProtoAccessor::template rebind<variable_element_type<V>>;
    return hana::type_c<Result>;
  }

  template <typename V>
  using accessor = typename decltype(get_accessor(hana::type_c<V>))::type;

  template <typename V> using value_type = typename accessor<V>::value_type;

  template <typename V> using element_type = typename accessor<V>::element_type;

  template <typename V> using pointer = typename accessor<V>::pointer;

  template <typename V>
  using const_pointer = typename accessor<std::add_const_t<V>>::pointer;

  template <typename V> using reference = typename accessor<V>::reference;

  template <typename V>
  using const_reference = typename accessor<std::add_const_t<V>>::reference;
};

///////////////////////////////////////////////////////////////////////////////
//                                                   variable_accessor_typedefs
//
// List case: We have multiple unions of variables at hand.

template <template <typename...> class List, typename ProtoAccessor,
          typename DefaultType, typename... Vs>
struct variable_accessor_typedefs<List<Vs...>, ProtoAccessor, DefaultType>
    : variable_find<List<Vs...>> {
  using base = variable_find<List<Vs...>>;
  using base::contains_collection;
  using base::find_collection_index;
  using base::get_variable_collections;

  template <typename V,
            typename = std::enable_if_t<
                find_collection_index(hana::type_c<V>) != hana::nothing>>
  static constexpr auto get_variable_accessor_typedefs(hana::basic_type<V>) {
    constexpr auto index = find_collection_index(hana::type_c<V>);
    constexpr auto unions = get_variable_collections();
    constexpr auto bucket = unions[*index];
    using W = typename decltype(bucket)::type;
    return variable_accessor_typedefs<W, ProtoAccessor, DefaultType>{};
  }

  template <typename V>
  using variable_accessor_typedefs_t =
      decltype(get_variable_accessor_typedefs(hana::type_c<V>));

  template <typename V>
  using accessor =
      typename variable_accessor_typedefs_t<V>::template accessor<V>;

  template <typename V>
  using value_type =
      typename variable_accessor_typedefs_t<V>::template value_type<V>;

  template <typename V>
  using element_type =
      typename variable_accessor_typedefs_t<V>::template element_type<V>;

  template <typename V>
  using pointer = typename variable_accessor_typedefs_t<V>::template pointer<V>;

  template <typename V>
  using const_pointer =
      typename variable_accessor_typedefs_t<V>::template const_pointer<V>;

  template <typename V>
  using reference =
      typename variable_accessor_typedefs_t<V>::template reference<V>;

  template <typename V>
  using const_reference =
      typename variable_accessor_typedefs_t<V>::template const_reference<V>;
};

///////////////////////////////////////////////////////////////////////////////
//                                             variable_accessor_implementation
//
// Base case: Only one variable union present.

template <typename VariableUnion, typename ProtoAccessor, typename DefaultType>
struct variable_accessor_implementation
    : variable_accessor_typedefs<VariableUnion, ProtoAccessor, DefaultType> {
  using base =
      variable_accessor_typedefs<VariableUnion, ProtoAccessor, DefaultType>;
  using base::contains_collection;
  using base::find_collection;
  using base::find_collection_index;
  using base::find_local_index;

  static constexpr auto get_value_storage() {
    using value_type_t = typename base::template value_type<VariableUnion>;
    constexpr std::size_t size = variable_traits<VariableUnion>::size();
    return hana::if_(hana::bool_c<size == 1>, hana::type_c<value_type_t>,
                     hana::type_c<std::array<value_type_t, size>>);
  }

  using base::pointer;

  using value_storage = typename decltype(get_value_storage())::type;
  using pointer_storage = typename base::template pointer<VariableUnion>;
  using const_pointer_storage =
      typename base::template const_pointer<VariableUnion>;

  template <typename V>
  using span = basic_span<typename base::template element_type<V>,
                          dynamic_extent, typename base::template accessor<V>>;

  template <typename>
  using const_span = basic_span<
      typename base::template element_type<std::add_const_t<VariableUnion>>,
      dynamic_extent,
      typename base::template accessor<std::add_const_t<VariableUnion>>>;

  /////////////////////////////////////////////////////////////////////////////
  //                      [variable_accessor_implementation.to_pointer_storage]

  template <typename Storage>
  static constexpr pointer_storage to_pointer_storage_impl(std::true_type,
                                                           Storage&& value) {
    return base::template accessor<VariableUnion>::to_pointer(value);
  }

  template <typename Storage>
  static constexpr pointer_storage to_pointer_storage_impl(std::false_type,
                                                           Storage&& value) {
    return base::template accessor<VariableUnion>::to_pointer(value[0]);
  }

  /// Transforms a value storage tuple into pointer storage tuple where each
  /// tuple entry points to the value tuple entry.
  static constexpr pointer_storage to_pointer_storage(value_storage& value) {
    return to_pointer_storage_impl(
        bool_c<variable_traits<VariableUnion>::size() == 1>, value);
  }

  /// Transforms a const value storage tuple into pointer to const storage tuple
  /// where each tuple entry points to the value tuple entry.
  static constexpr const_pointer_storage
  to_pointer_storage(const value_storage& value) {
    return to_pointer_storage_impl(
        bool_c<variable_traits<VariableUnion>::size() == 1>, value);
  }

  /////////////////////////////////////////////////////////////////////////////
  //                               [variable_accessor_implementation.make_span]

  template <typename V, typename = std::enable_if_t<
                            find_collection(hana::type_c<V>) != hana::nothing>>
  static constexpr span<V>
  make_span(hana::basic_type<V>, pointer_storage p,
            std::ptrdiff_t size_for_variable) noexcept {
    return span<V>(p,
                   variable_traits<VariableUnion>::size() * size_for_variable);
  }

  /////////////////////////////////////////////////////////////////////////////
  //                                  [variable_accessor_implementation.access]

  // Access a single Value

  template <typename V, typename S,
            typename = std::enable_if_t<contains_collection(hana::type_c<V>)>>
  static constexpr typename base::template reference<V>
  access_value(std::true_type, hana::basic_type<V>, S&& value) {
    auto p = base::template accessor<V>::to_pointer(value);
    return base::template accessor<V>::access(p, 1, 0);
  }

  template <typename V, typename S,
            typename = std::enable_if_t<find_collection(hana::type_c<V>) !=
                                        hana::nothing>>
  static constexpr typename base::template reference<V>
  access_value(std::false_type, hana::basic_type<V>, S&& value) {
    constexpr auto local_index = *find_local_index(hana::type_c<V>);
    auto p = base::template accessor<V>::to_pointer(value[local_index()]);
    return base::template accessor<V>::access(p, 1, 0);
  }

  template <typename V, typename = std::enable_if_t<
                            find_collection(hana::type_c<V>) != hana::nothing>>
  static constexpr typename base::template reference<V>
  access(hana::basic_type<V>, value_storage& value) {
    return access_value(bool_c<variable_traits<VariableUnion>::size() == 1>,
                        hana::type_c<V>, value);
  }

  template <typename V, typename = std::enable_if_t<
                            find_collection(hana::type_c<V>) != hana::nothing>>
  static constexpr typename base::template reference<std::add_const_t<V>>
  access(hana::basic_type<V>, const value_storage& value) {
    return access_value(bool_c<variable_traits<VariableUnion>::size() == 1>,
                        hana::type_c<std::add_const_t<V>>, value);
  }

  // Access a Span

  template <typename V, typename Span,
            typename = std::enable_if_t<is_span<Span>::value>,
            typename = std::enable_if_t<
                base::find_collection_index(hana::type_c<V>) != hana::nothing>>
  static constexpr decltype(auto) access_span_impl(hana::basic_type<V>,
                                                   const Span& span,
                                                   std::ptrdiff_t offset) {
    assert(variable_traits<VariableUnion>::size() <= span.size());
    assert(span.size() % variable_traits<VariableUnion>::size() == 0);
    std::ptrdiff_t size_for_variable =
        span.size() / variable_traits<VariableUnion>::size();
    constexpr auto variable_index = *find_local_index(hana::type_c<V>);
    return span[variable_index() * size_for_variable + offset];
  }

  /// Returns a reference to a variable value with specified storage, total
  /// extents and a given offset.
  template <typename V,
            typename = std::enable_if_t<
                find_collection_index(hana::type_c<V>) != hana::nothing>>
  static constexpr decltype(auto)
  access_const_span(hana::basic_type<V>,
                    nodeduce_t<const const_span<V>&> storage,
                    std::ptrdiff_t offset) {
    return access_span_impl(hana::type_c<V>, storage, offset);
  }

  /// Returns a reference to a variable value with specified storage, total
  /// extents and a given offset.
  template <typename V,
            typename = std::enable_if_t<
                find_collection_index(hana::type_c<V>) != hana::nothing>>
  static constexpr decltype(auto)
  access_span(hana::basic_type<V>, nodeduce_t<const span<V>&> storage,
              std::ptrdiff_t offset) {
    return access_span_impl(hana::type_c<V>, storage, offset);
  }

  template <typename V, typename S,
            typename = std::enable_if_t<find_collection(hana::type_c<V>) !=
                                        hana::nothing>,
            typename = std::enable_if_t<std::is_convertible<S, span<V>>::value>>
  static constexpr decltype(auto) access(hana::basic_type<V>, S&& storage,
                                         std::ptrdiff_t offset) {
    return access_span(hana::type_c<V>, storage, offset);
  }

  template <
      typename V, typename S,
      typename =
          std::enable_if_t<find_collection(hana::type_c<V>) != hana::nothing>,
      typename = std::enable_if_t<!std::is_convertible<S, span<V>>::value>,
      typename = std::enable_if_t<std::is_convertible<S, const_span<V>>::value>>
  static constexpr decltype(auto) access(hana::basic_type<V>, S&& storage,
                                         std::ptrdiff_t offset) {
    return access_const_span(hana::type_c<V>, storage, offset);
  }
};

///////////////////////////////////////////////////////////////////////////////
//                                             variable_accessor_implementation
//
// List case: Multiple unions.

template <template <typename...> class List, typename ProtoAccessor,
          typename DefaultType, typename... Vs>
struct variable_accessor_implementation<List<Vs...>, ProtoAccessor, DefaultType>
    : variable_accessor_typedefs<List<Vs...>, ProtoAccessor, DefaultType> {
  using base =
      variable_accessor_typedefs<List<Vs...>, ProtoAccessor, DefaultType>;
  using base::contains_collection;
  using base::find_collection;
  using base::find_collection_index;
  using base::get_variable_collections;

  template <typename V,
            typename = std::enable_if_t<
                find_collection_index(hana::type_c<V>) != hana::nothing>>
  static constexpr auto get_variable_accessor(hana::basic_type<V>) {
    constexpr auto index = find_collection_index(hana::type_c<V>);
    constexpr auto unions = get_variable_collections();
    constexpr auto bucket = unions[*index];
    using Union = typename decltype(bucket)::type;
    return variable_accessor_implementation<Union, ProtoAccessor,
                                            DefaultType>{};
  }

  template <typename V>
  using variable_accessor_t = decltype(get_variable_accessor(hana::type_c<V>));

  using value_storage =
      hana::tuple<typename variable_accessor_t<Vs>::value_storage...>;

  using pointer_storage =
      hana::tuple<typename variable_accessor_t<Vs>::pointer_storage...>;

  using const_pointer_storage =
      hana::tuple<typename variable_accessor_t<Vs>::const_pointer_storage...>;

  template <typename V>
  using span = typename variable_accessor_t<V>::template span<V>;

  template <typename V>
  using const_span = typename variable_accessor_t<V>::template const_span<V>;

  static constexpr pointer_storage to_pointer_storage(value_storage& storage) {
    constexpr auto unions = get_variable_collections();
    return hana::transform(unions, [&](auto vars) {
      constexpr auto union_ = decltype(vars){};
      constexpr auto union_index =
          *hana::index_if(unions, hana::equal.to(union_));
      return get_variable_accessor(union_).to_pointer_storage(
          storage[union_index]);
    });
  }

  static constexpr const_pointer_storage
  to_pointer_storage(const value_storage& storage) {
    constexpr auto unions = get_variable_collections();
    return hana::transform(unions, [&](auto vars) {
      constexpr auto add_const = hana::metafunction<std::add_const>;
      constexpr auto union_ = decltype(vars){};
      constexpr auto union_index =
          *hana::index_if(unions, hana::equal.to(union_));
      return get_variable_accessor(add_const(union_))
          .to_pointer_storage(storage[union_index]);
    });
  }

  template <typename V, typename = std::enable_if_t<
                            find_collection(hana::type_c<V>) != hana::nothing>>
  static constexpr span<V>
  make_span(hana::basic_type<V>, pointer_storage p,
            std::ptrdiff_t size_for_variable) noexcept {
    constexpr auto index = *find_collection_index(hana::type_c<V>);
    constexpr auto accessor = get_variable_accessor(hana::type_c<V>);
    return accessor.make_span(hana::type_c<V>, p[index], size_for_variable);
  }

  template <typename V, typename = std::enable_if_t<
                            find_collection(hana::type_c<V>) != hana::nothing>>
  static constexpr decltype(auto) access(hana::basic_type<V>,
                                         value_storage& storage) {
    constexpr auto index = *find_collection_index(hana::type_c<V>);
    constexpr auto accessor = get_variable_accessor(hana::type_c<V>);
    return accessor.access(hana::type_c<V>, storage[index]);
  }

  template <typename V, typename = std::enable_if_t<
                            find_collection(hana::type_c<V>) != hana::nothing>>
  static constexpr decltype(auto) access(hana::basic_type<V>,
                                         const value_storage& storage) {
    constexpr auto index = *find_collection_index(hana::type_c<V>);
    constexpr auto accessor = get_variable_accessor(hana::type_c<V>);
    return accessor.access(hana::type_c<V>, storage[index]);
  }

  template <typename V, typename = std::enable_if_t<
                            find_collection(hana::type_c<V>) != hana::nothing>>
  static constexpr decltype(auto) access(hana::basic_type<V>,
                                         nodeduce_t<const span<V>&> span,
                                         std::ptrdiff_t offset) {
    constexpr auto accessor = get_variable_accessor(hana::type_c<V>);
    return accessor.access(hana::type_c<V>, span, offset);
  }
};

} // namespace detail
} // namespace fub

#endif