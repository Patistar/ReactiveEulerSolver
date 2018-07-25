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

#ifndef FUB_GRID_QUANTITIES_REF_HPP
#define FUB_GRID_QUANTITIES_REF_HPP

#include "fub/quantities.hpp"

#include <boost/hana/for_each.hpp>

namespace fub {
namespace detail {
template <typename V> struct variables_remove_const; // {
//   using type = std::remove_const_t<V>;
// };
template <template <typename...> class L, typename... Vs>
struct variables_remove_const<L<Vs...>> {
  using type = L<std::remove_const_t<Vs>...>;
};

template <typename T>
using variables_remove_const_t = typename variables_remove_const<T>::type;
} // namespace detail

template <typename VariableList, typename Accessor>
class basic_quantities_ref : private Accessor {
public:
  using accessor = Accessor;

  using quantities_type =
      basic_quantities<detail::variables_remove_const_t<VariableList>,
                       typename Accessor::template rebind<
                           detail::variables_remove_const_t<VariableList>>>;

  template <typename V>
  using value_type = typename accessor::template value_type<V>;

  template <typename V>
  using reference = typename accessor::template reference<V>;

  using storage = typename accessor::pointer_storage;

  basic_quantities_ref() = default;

  // We want reference semantics, therefore non copyable.
  basic_quantities_ref(const basic_quantities_ref&) = delete;
  basic_quantities_ref& operator=(const basic_quantities_ref&) = delete;

  // Re-enable move assignable
  //   basic_quantities_ref(basic_quantities_ref&&) = default;
  //   basic_quantities_ref& operator=(basic_quantities_ref&&) = default;

  /// Constructs a quantities tuple from specified values.
  explicit constexpr basic_quantities_ref(const storage& s,
                                          const accessor& a = accessor())
      : Accessor(a), m_storage(s) {}

  basic_quantities_ref(quantities_type& q)
      : basic_quantities_ref(
            q.get_variable_accessor().to_pointer_storage(q.get_storage()),
            q.get_variable_accessor()) {}

  basic_quantities_ref(const quantities_type& q)
      : Accessor(q.get_variable_accessor()),
        m_storage(get_variable_accessor().to_pointer_storage(q.get_storage())) {
  }

  basic_quantities_ref& operator=(const quantities_type& q) {
    auto variables = q.get_variable_accessor().get_accessible_variables();
    hana::for_each(variables,
                   [&](auto variable) { get(variable) = q.get(variable); });
    return *this;
  }

  /// Constructs a quantities tuple from specified values.
  template <typename... Args>
  explicit constexpr basic_quantities_ref(accessor_arg_t, const accessor& a,
                                          Args&&... args)
      : basic_quantities_ref(storage(std::forward<Args>(args)...), a) {}

  /// Constructs a quantities tuple from specified values.
  //   template <typename Arg, typename... Args,
  //             typename = std::enable_if_t<
  //                 !std::is_same<remove_cvref_t<Arg>, storage>::value>>
  //   explicit constexpr basic_quantities_ref(Arg&& arg0, Args&&... args)
  //       : basic_quantities_ref(
  //             storage(std::forward<Arg>(arg0), std::forward<Args>(args)...))
  //             {}

  const accessor& get_variable_accessor() const noexcept { return *this; }

  const storage& get_storage() const noexcept { return m_storage; }

  template <typename V>
  reference<V> get(hana::basic_type<V> = hana::basic_type<V>()) const noexcept {
    const accessor& accessor = get_variable_accessor();
    return accessor.access(
        hana::type_c<V>, accessor.make_span(hana::type_c<V>, m_storage, 1), 0);
  }

  template <std::size_t I> decltype(auto) get() const noexcept {
    constexpr auto variables =
        get_variable_accessor().get_accessible_variables();
    using V = typename decltype(variables[size_c<I>])::type;
    return get<V>();
  }

  template <typename V> reference<V> operator[](const V&) const noexcept {
    return get<V>();
  }

  constexpr static std::size_t size() noexcept {
    return hana::size(Accessor::get_accessible_variables())();
  }

private:
  storage m_storage;
}; // namespace fub

template <typename... Variables>
using quantities_ref =
    basic_quantities_ref<types<Variables...>,
                         variable_accessor<types<Variables...>>>;

template <typename Variables, typename Abi = simd_abi::native<double>>
using simd_quantities_ref = basic_quantities_ref<
    Variables,
    variable_accessor<Variables, fub::accessor_simd_aligned<void, Abi>>>;

} // namespace fub

#endif