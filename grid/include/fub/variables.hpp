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

#include "fub/simd.hpp"
#include "fub/tuple.hpp"
#include "fub/type_traits.hpp"

#include <fmt/format.h>
#include <meta/meta.hpp>
#include <string>

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
}

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

  using pointer = std::conditional_t<std::is_const<V>::value, const value_type*,
                                     value_type*>;

  using const_pointer = const value_type*;

  using equal_to =
      typename detected_or<is_equal_fn<V>, equal_to_typedef, V>::type;

  static constexpr const char* name() noexcept { return V::name(); }

  /// @brief Returns the number of elements which are needed to cover the
  /// specified extents.
  template <typename Extents>
  static constexpr std::ptrdiff_t size(const Extents& extents) noexcept {
    return extents.size();
  }
};

template <typename V> struct variable_equal_fn {
  template <typename U>
  using invoke =
      disjunction<meta::invoke<typename variable_traits<V>::equal_to, U>,
                  meta::invoke<typename variable_traits<U>::equal_to, V>>;
};
// }}}

////////////////////////////////////////////////////////////////////////////////
//                                                    [meta.variable_find_index]

template <typename T, typename Tuple> struct variable_find_index;

template <typename T, template <typename...> class List, typename... Ts>
struct variable_find_index<T, List<Ts...>>
    : std::integral_constant<
          std::size_t,
          sizeof...(Ts) -
              meta::size<meta::find_if<meta::list<Ts...>,
                                       variable_equal_fn<T>>>::value> {};

template <typename T> struct variable_size;

template <template <typename...> class L, typename... Ts>
struct variable_size<L<Ts...>>
    : std::integral_constant<std::size_t, sizeof...(Ts)> {};

template <typename Var, typename Vars>
struct is_gettable
    : std::integral_constant<bool, (variable_find_index<Var, Vars>::value <
                                    variable_size<Vars>::value)> {};

////////////////////////////////////////////////////////////////////////////////
//                                                            [class.quantities]

template <typename... Vars> class quantities_ref;

/// @brief This class is a named tuple and provides access to data
/// assiciated via a variable tag.
template <typename... Vars> class quantities {
  static_assert(sizeof...(Vars) > 0, "There must be at least one variable.");
  std::tuple<typename variable_traits<Vars>::value_type...> m_values{};

public:
  quantities() noexcept = default;

  quantities(typename variable_traits<Vars>::value_type... values) noexcept
      : m_values{std::move(values)...} {}

  template <typename... OtherVars,
            typename = std::enable_if_t<conjunction<
                is_gettable<Vars, meta::list<OtherVars...>>...>::value>>
  quantities(const quantities<OtherVars...>& vars)
      : m_values{vars.template get<Vars>()...} {}

  template <typename... OtherVars,
            typename = std::enable_if_t<conjunction<
                is_gettable<Vars, meta::list<OtherVars...>>...>::value>>
  quantities(quantities_ref<OtherVars...> refs) noexcept
      : m_values{refs.template get<Vars>()...} {}

  template <typename Var> auto& get() noexcept {
    using Index = variable_find_index<Var, meta::list<Vars...>>;
    return std::get<Index::value>(m_values);
  }

  template <std::size_t I, typename = std::enable_if_t<(I < sizeof...(Vars))>>
  auto& get() noexcept {
    return std::get<I>(m_values);
  }

  template <typename Var> const auto& get() const noexcept {
    using Index = variable_find_index<Var, meta::list<Vars...>>;
    return std::get<Index::value>(m_values);
  }

  template <std::size_t I, typename = std::enable_if_t<(I < sizeof...(Vars))>>
  const auto& get() const noexcept {
    return std::get<I>(m_values);
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
  return std::tuple<Vars...>();
}

////////////////////////////////////////////////////////////////////////////////
//                                                         [class.quantities_ref]

/// @brief This class is a named tuple for references to data associated to
/// variables. It is useful to group some variables on a patch and modify it
/// there directly.
template <typename... Vars> class quantities_ref {
  std::tuple<typename variable_traits<Vars>::pointer...> m_pointers;

public:
  quantities_ref() = delete;
  quantities_ref(typename variable_traits<Vars>::pointer... pointers) noexcept
      : m_pointers{pointers...} {
    assert(ranges::tuple_foldl(m_pointers, true, std::logical_and<bool>()));
  }
  quantities_ref(const quantities_ref&) = default;

  quantities_ref& operator=(const quantities_ref& refs) {
    assign(refs);
    return *this;
  }

  template <typename... OtherVars>
  quantities_ref(quantities_ref<OtherVars...> data)
      : m_pointers{std::addressof(data.template get<Vars>())...} {}

  template <typename... OtherVars>
  quantities_ref& operator=(quantities_ref<OtherVars...> refs) noexcept {
    assign(refs);
    return *this;
  }

  template <typename... OtherVars>
  void assign(quantities_ref<OtherVars...> refs) noexcept {
    fub::for_each_tuple_element(
        [&, this](auto var) {
          using T = decltype(var);
          this->get<T>() = refs.template get<T>();
        },
        Vars{}...);
  }

  /// Assign data from a given variable set
  template <typename... OtherVars>
  quantities_ref& operator=(const quantities<OtherVars...>& vars) noexcept {
    fub::for_each_tuple_element(
        [&, this](auto var) {
          using T = decltype(var);
          this->get<T>() = vars.template get<T>();
        },
        Vars{}...);
    return *this;
  }

  /// @brief Returns a reference to the data associated with variable Var.
  template <typename Var> auto& get() const noexcept {
    using I = variable_find_index<std::remove_const_t<Var>,
                                  meta::list<std::remove_const_t<Vars>...>>;
    return *std::get<I::value>(m_pointers);
  }
  /// @brief Returns a reference to the data associated with variable Var.
  template <std::size_t I> auto& get() const noexcept {
    return *std::get<I>(m_pointers);
  }

  /// @brief Returns a reference to the data associated with variable Var.
  template <typename Var> auto& operator[](Var) const noexcept {
    return get<Var>();
  }
};

template <typename... Vars>
auto get_variables(const quantities_ref<Vars...>&) noexcept {
  return std::tuple<Vars...>();
}

////////////////////////////////////////////////////////////////////////////////
//                                                           [meta.flux_variable]

template <typename Var> struct flux_variable;

template <typename V> struct is_flux_variable : std::false_type {};
template <typename V>
struct is_flux_variable<flux_variable<V>> : std::true_type {};

/// @brief Applies the flux_variable function on or many quantities.
template <typename Variable> struct add_flux {
  using type = flux_variable<Variable>;
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

/// @brief This meta function is used to transform a concrete variable to a flux
/// variable.
///
/// This saves boilerplate code as it automatically creates a distinct variable
/// by prepending the string "Flux " for the base variable name.
template <typename Var> struct flux_variable {
  static const std::string name_;
  using value_type = typename variable_traits<Var>::value_type;
  static const char* name() noexcept { return name_.c_str(); }

  struct equal_to {
    template <typename U>
    using invoke = disjunction<
        std::is_same<flux_variable, std::remove_cv_t<U>>,
        conjunction<is_flux_variable<U>,
                    meta::invoke<typename variable_traits<Var>::equal_to,
                                 remove_flux_t<U>>>>;
  };
};

template <typename Var>
const std::string flux_variable<Var>::name_{
    fmt::format("Flux {}", Var::name())};

////////////////////////////////////////////////////////////////////////////////
//                                                               [meta.MakeFlux]
//                                                                   [meta.Flux]

template <typename Var> flux_variable<Var> flux(Var) { return {}; }
template <typename Var> Var unflux(flux_variable<Var>) { return {}; }

////////////////////////////////////////////////////////////////////////////////
//                                                           [meta.simd_wrapper]

template <typename Variable,
          typename Abi = simd_abi::native<typename Variable::value_type>>
struct simd_wrapper {
  static constexpr decltype(auto) name() noexcept { return Variable::name(); }
  using value_type = simd<typename variable_traits<Variable>::value_type, Abi>;
  struct equal_to {
    template <typename U>
    using invoke = disjunction<
        std::is_same<simd_wrapper, U>,
        meta::invoke<typename variable_traits<Variable>::equal_to, U>>;
  };
};

template <typename V, typename Abi = simd_abi::native<double>> struct add_simd {
  using type = simd_wrapper<V, Abi>;
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
} // namespace fub

/// @brief These specialisations make stuctured binding work with quantities.
///
/// Example:
///    quantities<Density, Velocity> var{1., 2.}
///    auto [rho, u] = var;
///    assert(rho == 1.);
///    assert(u  == 2.);
///
/// @{
namespace std {
template <typename... Vars>
class tuple_size<::fub::quantities<Vars...>>
    : public integral_constant<std::size_t, sizeof...(Vars)> {};

template <std::size_t I, typename... Vars>
class tuple_element<I, ::fub::quantities<Vars...>>
    : public tuple_element<
          I, std::tuple<typename ::fub::variable_traits<Vars>::value_type...>> {
};

template <typename... Vars>
class tuple_size<::fub::quantities_ref<Vars...>>
    : public integral_constant<std::size_t, sizeof...(Vars)> {};

template <std::size_t I, typename... Vars>
class tuple_element<I, ::fub::quantities_ref<Vars...>>
    : public tuple_element<
          I, std::tuple<typename ::fub::variable_traits<Vars>::value_type...>> {
};
} // namespace std
/// @}

#endif
