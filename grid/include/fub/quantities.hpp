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

#ifndef FUB_GRID_QUANTITIES_HPP
#define FUB_GRID_QUANTITIES_HPP

#include "fub/accessor_simd.hpp"
#include "fub/variable_accessor.hpp"

namespace fub {
struct accessor_arg_t {};
static constexpr accessor_arg_t accessor_arg{};

template <typename VariableList, typename Accessor>
class basic_quantities : private Accessor {
public:
  static_assert(
      hana::not_(hana::any_of(hana::to_tuple(VariableList{}),
                              hana::trait<std::is_const>)),
      "basic_quantities<...> is a value type. No variable should be const");

  using accessor = Accessor;
  using variable_list = VariableList;

  template <typename V>
  using value_type = typename accessor::template value_type<V>;

  template <typename V>
  using reference = typename accessor::template reference<V>;

  template <typename V>
  using const_reference =
      typename accessor::template reference<std::add_const_t<V>>;

  using storage = typename accessor::value_storage;

  basic_quantities() = default;
  basic_quantities(const basic_quantities&) = default;
  basic_quantities(basic_quantities&&) = default;
  basic_quantities& operator=(const basic_quantities&) = default;
  basic_quantities& operator=(basic_quantities&&) = default;

  /// Constructs a quantities tuple from specified values.
  explicit constexpr basic_quantities(const storage& s,
                                      const accessor& a = accessor())
      : Accessor(a), m_storage(s) {}

  /// Constructs a quantities tuple from specified values.
  template <typename... Args>
  explicit constexpr basic_quantities(accessor_arg_t, const accessor& a,
                                      Args&&... args)
      : basic_quantities(storage(std::forward<Args>(args)...), a) {}

  /// Constructs a quantities tuple from specified values.
  template <typename Arg, typename... Args,
            typename = std::enable_if_t<
                !std::is_same<remove_cvref_t<Arg>, storage>::value>>
  explicit constexpr basic_quantities(Arg&& arg0, Args&&... args)
      : basic_quantities(
            storage(std::forward<Arg>(arg0), std::forward<Args>(args)...)) {}

  const accessor& get_variable_accessor() const noexcept { return *this; }

  const storage& get_storage() const noexcept { return m_storage; }
  storage& get_storage() noexcept { return m_storage; }

  template <std::size_t I> decltype(auto) get() noexcept {
    constexpr auto variables =
        get_variable_accessor().get_accessible_variables();
    return get_variable_accessor().access(variables[size_c<I>], m_storage);
  }

  template <std::size_t I> decltype(auto) get() const noexcept {
    constexpr auto variables =
        get_variable_accessor().get_accessible_variables();
    return get_variable_accessor().access(variables[size_c<I>], m_storage);
  }

  template <typename V>
  reference<V> get(hana::basic_type<V> = hana::basic_type<V>()) noexcept {
    return get_variable_accessor().access(hana::type_c<V>, m_storage);
  }

  template <typename V>
  const_reference<V> get(hana::basic_type<V> = hana::basic_type<V>()) const
      noexcept {
    return get_variable_accessor().access(hana::type_c<V>, m_storage);
  }

  template <typename V> reference<V> operator[](const V&) noexcept {
    return get<V>();
  }

  template <typename V> const_reference<V> operator[](const V&) const noexcept {
    return get<V>();
  }

  constexpr static std::size_t size() noexcept {
    return hana::size(Accessor::get_accessible_variables())();
  }

  friend bool operator==(const basic_quantities& lhs,
                         const basic_quantities& rhs) noexcept {
    return lhs.m_storage == rhs.m_storage;
  }

  friend bool operator!=(const basic_quantities& lhs,
                         const basic_quantities& rhs) noexcept {
    return lhs.m_storage != rhs.m_storage;
  }

private:
  storage m_storage;
}; // namespace fub

template <typename... Variables>
using quantities = basic_quantities<types<Variables...>,
                                    variable_accessor<types<Variables...>>>;

template <typename Variables, typename Abi = simd_abi::native<double>>
using simd_quantities = basic_quantities<
    Variables,
    variable_accessor<Variables, fub::accessor_simd_aligned<void, Abi>>>;

} // namespace fub

#endif