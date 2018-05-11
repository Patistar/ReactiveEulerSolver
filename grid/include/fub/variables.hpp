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

#include "fub/extents.hpp"
#include "fub/mdspan.hpp"
#include "fub/simd.hpp"
#include "fub/tuple.hpp"
#include "fub/type_traits.hpp"

#include <fmt/format.h>
#include <meta/meta.hpp>
#include <string>
#include <vector>

namespace fub {

////////////////////////////////////////////////////////////////////////////////
//                                                       [class.variable_traits]
// {{{

template <typename V> using equal_to_typedef = typename V::equal_to;

template <typename V> struct is_equal_fn {
  template <typename U>
  using invoke = std::is_same<std::remove_cv_t<V>, std::remove_cv_t<U>>;
};

namespace detail {
template <typename V> using value_type_t = typename V::value_type;
template <typename V> using pointer_t = typename V::pointer;
template <typename V> using const_pointer_t = typename V::const_pointer;
template <typename V> using variables_tuple_t = typename V::variables_tuple;

template <typename V, typename A>
using storage_type_t = typename V::template storage_type<A>;

template <typename V, typename S>
using variable_get_pointer_memfn_t =
    decltype(V::get_pointer(std::declval<S>(), 0));

template <typename V>
using variable_size_memfn_t = decltype(V::size(extents<1>()));

template <typename V, typename W, typename T>
using variable_access_memfn_t = decltype(V::access(W(), std::declval<T>()));

template <typename V, typename W, typename P, typename E>
using variable_view_memfn_t = decltype(V::view(W(), P(), std::declval<E>()));

} // namespace detail

/// @brief The `variable_traits` class template provides uniform access to
/// various properties of variables.
///
/// The class templates `quantities`, `patch` and `patch_view` use this helper
/// class to access these properties.
template <typename V> struct variable_traits {
  /// @brief The variable type
  using variable_type = V;

  /// @brief The underlying value type, for example `double`, `float` or `int`.
  using value_type =
      typename detected_or<double, detail::value_type_t, V>::type;

  template <typename Allocator>
  using storage_type =
      typename detected_or<std::vector<value_type, Allocator>,
                           detail::storage_type_t, V, Allocator>::type;

  using pointer = typename detected_or<value_type*, detail::pointer_t, V>::type;

  using const_pointer =
      typename detected_or<const value_type*, detail::const_pointer_t, V>::type;

  using equal_to =
      typename detected_or<is_equal_fn<V>, equal_to_typedef, V>::type;

  static constexpr const char* name() noexcept { return V::name(); }

  /// @brief Returns the number of elements which are needed to cover the
  /// specified extents.
  template <typename Extents>
  static constexpr std::ptrdiff_t size(const Extents& extents) noexcept {
    if constexpr (is_detected<detail::variable_size_memfn_t, V>::value) {
      return V::size(extents);
    } else {
      return extents.size();
    }
  }

  template <typename W, typename T>
  static constexpr decltype(auto) access(W w, T&& x) noexcept {
    if constexpr (is_detected<detail::variable_access_memfn_t, V, W,
                              T>::value) {
      return V::access(w, x);
    } else if constexpr (std::is_pointer<std::remove_reference_t<T>>::value) {
      return (*x);
    } else {
      return std::forward<T>(x);
    }
  }

  template <typename W, typename E>
  static constexpr decltype(auto) view(W w, pointer ptr,
                                       const E& extents) noexcept {
    if constexpr (is_detected<detail::variable_view_memfn_t, V, W, pointer,
                              const E&>::value) {
      return V::view(w, ptr, extents);
    } else {
      return mdspan<value_type, E>{ptr, extents};
    }
  }

  template <typename W, typename E>
  static constexpr decltype(auto) view(W w, const_pointer ptr,
                                       const E& extents) noexcept {
    if constexpr (is_detected<detail::variable_view_memfn_t, V, W,
                              const_pointer, const E&>::value) {
      return V::view(w, ptr, extents);
    } else {
      return mdspan<const value_type, E>{ptr, extents};
    }
  }

  template <typename S>
  static constexpr auto get_pointer(S&& storage, index i = 0) {
    if constexpr (is_detected<detail::variable_get_pointer_memfn_t, V,
                              S>::value) {
      return V::get_pointer(std::forward<S>(storage), i);
    } else {
      return storage.data() + i;
    }
  }
};

template <typename V> struct variable_equal_fn {
  template <typename U>
  using invoke =
      disjunction<std::is_same<U, V>,
                  meta::invoke<typename variable_traits<V>::equal_to, U>,
                  meta::invoke<typename variable_traits<U>::equal_to, V>>;
};
// }}}

////////////////////////////////////////////////////////////////////////////////
//                                                    [meta.variable_find_index]
// {{{
template <typename T, typename Tuple> struct variable_find_index;

template <typename T, template <typename...> class List, typename... Ts>
struct variable_find_index<T, List<Ts...>>
    : std::integral_constant<
          std::size_t,
          sizeof...(Ts) -
              meta::size<meta::find_if<
                  meta::list<std::remove_const_t<Ts>...>,
                  variable_equal_fn<std::remove_const_t<T>>>>::value> {};

template <typename T> struct variable_size;

template <template <typename...> class L, typename... Ts>
struct variable_size<L<Ts...>>
    : std::integral_constant<std::size_t, sizeof...(Ts)> {};

template <typename Var, typename Vars>
struct is_gettable
    : std::integral_constant<bool, (variable_find_index<Var, Vars>::value <
                                    variable_size<Vars>::value)> {};
// }}}

////////////////////////////////////////////////////////////////////////////////
//                                                  [function.flatten_variables]
// {{{
template <typename V> auto flatten_variables(V) {
  if constexpr (is_detected<detail::variables_tuple_t, V>::value) {
    return typename V::variables_tuple{};

  } else {
    return std::tuple<V>{};
  }
}

template <typename Var> struct flux_variable;

template <typename V, typename... Vars>
auto flatten_variables(flux_variable<V>, Vars... vars) {
  return flatten_variables(V{}, vars...);
}

template <typename V, typename... Vars> auto flatten_variables(V, Vars...) {
  if constexpr (is_detected<detail::variables_tuple_t, V>::value) {
    return std::tuple_cat(typename V::variables_tuple{},
                          flatten_variables(Vars{}...));

  } else {
    return std::tuple_cat(std::tuple<V>{}, flatten_variables(Vars{}...));
  }
}

template <typename... Vars>
using flatten_variables_t =
    decltype(flatten_variables(std::declval<Vars>()...));
// }}}

////////////////////////////////////////////////////////////////////////////////
//                                                        [meta.vector_variable]
// {{{
template <typename V, typename... Vars> struct vector_variable {
  using variables_tuple = std::tuple<V, Vars...>;
  using scalar_type = typename variable_traits<V>::value_type;
  static constexpr int rank = 1 + sizeof...(Vars);

  struct equal_to {
    template <typename U>
    using invoke = disjunction<
        std::is_same<vector_variable, U>,
        meta::invoke<typename variable_traits<V>::equal_to, U>,
        meta::invoke<typename variable_traits<Vars>::equal_to, U>...>;
  };

  static_assert(
      conjunction<std::is_same<
          scalar_type, typename variable_traits<Vars>::value_type>...>::value);
  using value_type = std::array<scalar_type, rank>;

  struct pointer {
    scalar_type* native;
    index stride;

    operator bool() const noexcept { return native; }

    friend pointer operator+(pointer p, index offset) noexcept {
      return {p.native + offset, p.stride};
    }
  };

  template <typename Alloc>
  static pointer get_pointer(std::vector<scalar_type, Alloc>& store,
                             index i = 0) noexcept {
    const index size = store.size() / rank;
    return {store.data() + i, size};
  }

  template <typename Extents>
  static pointer get_pointer(const mdspan<scalar_type, Extents>& view,
                             index i = 0) noexcept {
    return {view.data() + i, view.size()};
  }

  struct const_pointer {
    const scalar_type* native;
    index stride;

    operator bool() const noexcept { return native; }

    friend const_pointer operator+(const_pointer cp, index offset) noexcept {
      return {cp.native + offset, cp.stride};
    }
  };

  template <typename Alloc>
  static const_pointer get_pointer(const std::vector<scalar_type, Alloc>& store,
                                   index i = 0) noexcept {
    const index size = store.size() / rank;
    return {store.data() + i, size};
  }

  template <typename Extents>
  static const_pointer
  get_pointer(const mdspan<const scalar_type, Extents>& view,
              index i = 0) noexcept {
    return {view.data() + i, view.size()};
  }

  template <typename W>
  static constexpr decltype(auto) access(W, pointer ptr) noexcept {
    using I = variable_find_index<W, meta::list<V, Vars...>>;
    static_assert(I::value < rank);
    return *(ptr.native + I::value * ptr.stride);
  }

  template <typename W>
  static constexpr decltype(auto) access(W, const_pointer ptr) noexcept {
    using I = variable_find_index<W, meta::list<V, Vars...>>;
    static_assert(I::value < rank);
    return *(ptr.native + I::value * ptr.stride);
  }

  template <typename W>
  static constexpr decltype(auto) access(W, value_type& array) noexcept {
    using I = variable_find_index<W, meta::list<V, Vars...>>;
    static_assert(I::value < rank);
    return array[I::value];
  }

  template <typename W>
  static constexpr decltype(auto) access(W, const value_type& array) noexcept {
    using I = variable_find_index<W, meta::list<V, Vars...>>;
    static_assert(I::value < rank);
    return array[I::value];
  }

  // template <typename Allocator>
  // using storage_type = std::vector<scalar_type, Allocator>;

  /// @brief Returns the number of elements which are needed to cover the
  /// specified extents.
  template <typename Extents>
  static constexpr index size(const Extents& extents) noexcept {
    return rank * variable_traits<V>::size(extents);
  }

  template <typename W, typename Extents>
  static mdspan<scalar_type, Extents> view(W, pointer p,
                                           const Extents& extents) noexcept {
    using I = variable_find_index<W, meta::list<V, Vars...>>;
    static_assert(I::value < rank);
    return mdspan<scalar_type, Extents>{p.native + I::value * p.stride,
                                        extents};
  }

  template <typename W, typename Extents>
  static mdspan<const scalar_type, Extents>
  view(W, const_pointer p, const Extents& extents) noexcept {
    using I = variable_find_index<W, meta::list<V, Vars...>>;
    static_assert(I::value < rank);
    return mdspan<const scalar_type, Extents>{p.native + I::value * p.stride,
                                              extents};
  }
};

template <typename T>
struct is_vector_variable : is_detected<detail::variables_tuple_t, T> {};
// }}}

////////////////////////////////////////////////////////////////////////////////
//                                                            [class.quantities]
// {{{

template <typename... Vars> class quantities_ref;

/// @brief This class is a named tuple and provides access to data
/// assiciated via a variable tag.
template <typename... Vars> class quantities {
  static_assert(sizeof...(Vars) > 0, "There must be at least one variable.");
  std::tuple<typename variable_traits<Vars>::value_type...> m_values{};

  template <typename Var, typename T> auto construct(Var, T&& source) {
    if constexpr (is_vector_variable<Var>::value) {
      return fub::apply(
          [&](auto... vars) {
            return typename variable_traits<Var>::value_type{source[vars]...};
          },
          typename Var::variables_tuple{});
    } else {
      return source.template get<Var>();
    }
  }

public:
  quantities() noexcept = default;

  quantities(typename variable_traits<Vars>::value_type... values) noexcept
      : m_values{std::move(values)...} {}

  template <typename... OtherVars>
  quantities(const quantities<OtherVars...>& vars)
      : m_values{construct(Vars{}, vars)...} {}

  template <typename... OtherVars>
  quantities(quantities_ref<OtherVars...> refs)
      : m_values{construct(Vars{}, refs)...} {}

  template <typename Var> auto& get() noexcept {
    using I = variable_find_index<Var, meta::list<Vars...>>;
    using Holder = std::tuple_element_t<I::value, std::tuple<Vars...>>;
    return variable_traits<Holder>::access(Var(), std::get<I::value>(m_values));
  }

  template <typename Var> const auto& get() const noexcept {
    using I = variable_find_index<Var, meta::list<Vars...>>;
    using Holder = std::tuple_element_t<I::value, std::tuple<Vars...>>;
    return variable_traits<Holder>::access(Var(), std::get<I::value>(m_values));
  }

  template <typename Var> auto& operator[](Var) noexcept { return get<Var>(); }

  template <typename Var> const auto& operator[](Var) const noexcept {
    return get<Var>();
  }
};

template <> class quantities<> {};

template <typename List> struct as_quantities : list_cast<quantities, List> {};
template <typename List>
using as_quantities_t = typename as_quantities<List>::type;

template <typename... Vars>
auto get_variables(
    const quantities<Vars...>& = quantities<Vars...>()) noexcept {
  return std::tuple<Vars...>{};
}
// }}}

////////////////////////////////////////////////////////////////////////////////
//                                                        [class.quantities_ref]
// {{{

/// @brief This class is a named tuple for references to data associated to
/// variables. It is useful to group some variables on a patch and modify it
/// there directly.
template <typename... Vars> class quantities_ref {
  template <typename V>
  using pointer = std::conditional_t<
      std::is_const<V>::value,
      typename variable_traits<std::remove_const_t<V>>::const_pointer,
      typename variable_traits<V>::pointer>;
  std::tuple<pointer<Vars>...> m_pointers;

  template <typename... Wars> friend class quantities_ref;

public:
  quantities_ref() = delete;
  quantities_ref(pointer<Vars>... pointers) noexcept : m_pointers{pointers...} {
    assert(true && ... && pointers);
  }
  quantities_ref(const quantities_ref&) = default;

  quantities_ref& operator=(const quantities_ref& refs) {
    assign(refs);
    return *this;
  }

  template <typename... OtherVars>
  quantities_ref(quantities_ref<OtherVars...> other)
      : m_pointers{other.m_pointers} {}

  template <typename... OtherVars>
  quantities_ref& operator=(quantities_ref<OtherVars...> refs) noexcept {
    assign(refs);
    return *this;
  }

  template <typename... OtherVars>
  void assign(const quantities_ref<OtherVars...>& refs) noexcept {
    for_each_tuple_element(
        [&, this](auto var) {
          using T = decltype(var);
          this->get<T>() = refs.template get<T>();
        },
        flatten_variables(Vars{}...));
  }

  /// Assign data from a given variable set
  template <typename... OtherVars>
  quantities_ref& operator=(const quantities<OtherVars...>& vars) noexcept {
    fub::for_each_tuple_element(
        [&, this](auto var) {
          using T = decltype(var);
          this->get<T>() = vars.template get<T>();
        },
        flatten_variables(Vars{}...));
    return *this;
  }

  /// @brief Returns a reference to the data associated with variable Var.
  template <typename Var> auto& get() const noexcept {
    using I = variable_find_index<Var, meta::list<Vars...>>;
    using Holder = std::tuple_element_t<I::value, std::tuple<Vars...>>;
    using std::get;
    return variable_traits<Holder>::access(Var(), get<I::value>(m_pointers));
  }

  /// @brief Returns a reference to the data associated with variable Var.
  template <typename Var> auto& operator[](Var) const noexcept {
    return get<Var>();
  }
};

template <typename... Vars>
auto get_variables(const quantities_ref<Vars...>&) noexcept {
  return flatten_variables(Vars{}...);
}
// }}}

////////////////////////////////////////////////////////////////////////////////
//                                                          [meta.flux_variable]
// {{{
template <typename V> struct is_flux_variable : std::false_type {};
template <typename V>
struct is_flux_variable<flux_variable<V>> : std::true_type {};

template <typename Variable> struct add_flux;

namespace detail {
struct fluxify {
  template <typename W> using invoke = typename add_flux<W>::type;
};

template <typename V, bool> struct add_flux_;

template <typename V> struct add_flux_<V, false> {
  using type = flux_variable<V>;
};

template <typename V> struct add_flux_<V, true> {
  struct type
      : list_cast_t<vector_variable,
                    meta::transform<
                        list_cast_t<meta::list, typename V::variables_tuple>,
                        fluxify>> {};
};
} // namespace detail

/// @brief Applies the flux_variable function on or many quantities.
template <typename Variable> struct add_flux {
  using type =
      typename detail::add_flux_<Variable,
                                 is_vector_variable<Variable>::value>::type;
};
template <typename V> using add_flux_t = typename add_flux<V>::type;

/// @brief Specialisation for a set of variables.
template <typename... Vars> struct add_flux<std::tuple<Vars...>> {
  using type = std::tuple<add_flux_t<Vars>...>;
};
template <typename... Vars> struct add_flux<quantities<Vars...>> {
  using type = quantities<add_flux_t<Vars>...>;
};

template <typename Var> struct remove_flux { using type = Var; };
template <typename Var> using remove_flux_t = typename remove_flux<Var>::type;
template <typename Var> struct remove_flux<flux_variable<Var>> {
  using type = Var;
};
template <typename... Vars> struct remove_flux<std::tuple<Vars...>> {
  using type = std::tuple<remove_flux_t<Vars>...>;
};
template <typename... Vars> struct remove_flux<quantities<Vars...>> {
  using type = quantities<remove_flux_t<Vars>...>;
};

/// @brief This meta function is used to transform a concrete variable to a
/// flux variable.
///
/// This saves boilerplate code as it automatically creates a distinct
/// variable by prepending the string "Flux " for the base variable name.
template <typename Var> struct flux_variable {
  static const std::string name_;
  using value_type = typename variable_traits<Var>::value_type;
  static const char* name() noexcept { return name_.c_str(); }

  struct equal_to {
    template <typename U>
    using invoke =
        disjunction<std::is_same<flux_variable, std::remove_cv_t<U>>,
                    meta::invoke<typename variable_traits<Var>::equal_to,
                                 remove_flux_t<U>>>;
  };

  template <typename W, typename P> decltype(auto) access(W w, P&& ptr) {
    return variable_traits<Var>::access(w, std::forward<P>(ptr));
  }
};

template <typename Var>
const std::string flux_variable<Var>::name_{
    fmt::format("Flux {}", Var::name())};

template <typename Var> flux_variable<Var> flux(Var) { return {}; }
template <typename Var> remove_flux_t<Var> unflux(Var) {
  return remove_flux_t<Var>{};
}
// }}}

////////////////////////////////////////////////////////////////////////////////
//                                                           [meta.simd_wrapper]
// {{{
template <typename Variable,
          typename Abi =
              simd_abi::native<typename variable_traits<Variable>::value_type>>
struct simd_wrapper {
  static constexpr decltype(auto) name() noexcept { return Variable::name(); }
  using value_type = simd<typename variable_traits<Variable>::value_type, Abi>;
  struct equal_to {
    template <typename U>
    using invoke = disjunction<
        std::is_same<simd_wrapper, U>,
        meta::invoke<typename variable_traits<Variable>::equal_to, U>>;
  };
  template <typename W, typename T> decltype(auto) access(W w, T&& x) {
    return variable_traits<Variable>::access(w, std::forward<T>(x));
  }
};

template <typename V, typename Abi = simd_abi::native<double>> struct add_simd;

namespace detail {
template <typename Abi> struct simdify {
  template <typename W> using invoke = typename add_simd<W, Abi>::type;
};

template <typename V, typename Abi, bool> struct add_simd_;

template <typename V, typename Abi> struct add_simd_<V, Abi, false> {
  using type = simd_wrapper<V, Abi>;
};

template <typename V, typename Abi> struct add_simd_<V, Abi, true> {
  struct type
      : list_cast_t<vector_variable,
                    meta::transform<
                        list_cast_t<meta::list, typename V::variables_tuple>,
                        simdify<Abi>>> {};
};
} // namespace detail

template <typename V, typename Abi> struct add_simd {
  using type =
      typename detail::add_simd_<V, Abi, is_vector_variable<V>::value>::type;
};

template <typename V, typename Abi = simd_abi::native<double>>
using add_simd_t = typename add_simd<V, Abi>::type;

template <typename V>
using add_scalar_t = typename add_simd<V, simd_abi::scalar>::type;
template <typename V, typename Abi> struct add_simd<simd_wrapper<V>, Abi> {
  using type = simd_wrapper<V, Abi>;
};
template <typename Abi, typename... Vs>
struct add_simd<quantities<Vs...>, Abi> {
  using type = quantities<add_simd_t<Vs, Abi>...>;
};
template <typename V, typename Abi> struct add_simd<flux_variable<V>, Abi> {
  using type = flux_variable<add_simd_t<V, Abi>>;
};
template <typename V, typename Abi = simd_abi::native<double>>
add_simd_t<V, Abi> as_simd(V) {
  return {};
}
template <typename V> V unsimd(simd_wrapper<V>) { return {}; }

template <typename V>
using physical_dimensions_t =
    decltype(std::declval<const V&>().physical_dimensions());

template <typename V>
struct variable_is_dimensional : is_detected<physical_dimensions_t, V> {};
// }}}
} // namespace fub
/// @}

#endif
