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

#ifndef FUB_GRID_VARIABLE_HPP
#define FUB_GRID_VARIABLE_HPP

#include "fub/extents.hpp"
#include "fub/functional.hpp"
#include "fub/hana.hpp"
#include "fub/mdspan.hpp"
#include "fub/optional.hpp"
#include "fub/type_traits.hpp"
#include "fub/variant.hpp"

#include <boost/hana.hpp>

namespace fub {
inline namespace v1 {
namespace hana = boost::hana;

namespace detail {
template <typename T> using static_size_fn = decltype(T::static_size());
template <typename T> using size_memfn = decltype(std::declval<T>().size());
} // namespace detail

/// This is a trait type which checks if the `T` fulfills the \b Variable
/// requirements.
///
/// \see vector_variable
template <typename T>
struct is_variable
    : conjunction<is_regular<T>, is_detected<detail::static_size_fn, T>,
                  is_detected<detail::size_memfn, const T&>> {};

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
    return extents<Rank>::static_extent(0);
  }

  /// Returns the run-time Rank of the variable.
  ///
  /// \throws Nothing.
  constexpr std::ptrdiff_t size() const noexcept { return extents<Rank>::extent(0); }
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

namespace detail {
// This trait checks if `N` is a valid size
template <std::ptrdiff_t N, std::ptrdiff_t Size>
struct is_both_dynamic_or_less_than
    : std::integral_constant<
          bool, ((N == dynamic_extent && Size == dynamic_extent) || N < Size)> {
};
} // namespace detail

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

  static_assert(
      detail::is_both_dynamic_or_less_than<N, V::static_size()>::value,
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

} // namespace v1
} // namespace fub

#endif