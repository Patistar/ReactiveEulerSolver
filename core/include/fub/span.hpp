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

#ifndef FUB_CORE_SPAN_HPP
#define FUB_CORE_SPAN_HPP

#include "fub/array.hpp"
#include "fub/type_traits.hpp"
#include <fmt/format.h>
#include <vector>

namespace fub {

template <typename T, index N = dyn> class span {
public:
  // Types

  static_assert(N > 0);
  using element_type = T;
  using value_type = std::remove_const_t<T>;
  using pointer = T*;
  using reference = T&;
  using iterator = pointer;
  using const_iterator = const T*;

  // Constructors

  constexpr span() = default;

  template <std::size_t M, typename = std::enable_if_t<(N <= M)>>
  constexpr span(element_type (&a)[M]) noexcept : m_pointer{a.data()} {}

  span(pointer ptr, index size) : m_pointer{ptr} {
    if (size < N) {
      std::string what =
          fmt::format("span::span: Attempted to create a span of length '{}' "
                      "with a specified vector of length '{}'.",
                      N, size);
      throw std::out_of_range{what};
    }
  }

  template <
      typename S, index M,
      typename = std::enable_if_t<std::is_convertible<S*, pointer>::value>,
      typename = std::enable_if_t<(N <= M)>>
  constexpr span(const span<S, M>& other) noexcept : m_pointer{other.data()} {}

  template <
      typename S, index M,
      typename = std::enable_if_t<std::is_convertible<S*, pointer>::value>,
      typename = std::enable_if_t<(N <= M)>>
  constexpr span operator=(const span<S, M>& other) noexcept {
    m_pointer = other.data();
    return *this;
  }

  template <std::size_t M, typename S = element_type,
            typename = std::enable_if_t<std::is_const<S>::value && (N <= M)>>
  constexpr span(const array<std::remove_const_t<element_type>, M>& a) noexcept
      : m_pointer{a.data()} {}

  template <std::size_t M, typename = std::enable_if_t<(N <= M)>>
  constexpr span(array<std::remove_const_t<element_type>, M>& a) noexcept
      : m_pointer{a.data()} {}

  template <typename Alloc>
  span(std::vector<std::remove_const_t<element_type>, Alloc>& v)
      : m_pointer{v.data()} {
    if (v.size() < N) {
      std::string what =
          fmt::format("span::span: Attempted to create a span of length '{}' "
                      "with a specified vector of length '{}'.",
                      N, v.size());
      throw std::out_of_range{what};
    }
  }

  template <typename Alloc, typename S = element_type,
            typename = std::enable_if_t<std::is_const<S>::value>>
  span(const std::vector<std::remove_const_t<element_type>, Alloc>& v)
      : m_pointer{v.data()} {
    if (v.size() < N) {
      std::string what =
          fmt::format("span::span: Attempted to create a span of length '{}' "
                      "with a specified vector of length '{}'.",
                      N, v.size());
      throw std::out_of_range{what};
    }
  }

  // Class Member Access

  constexpr pointer data() const noexcept { return m_pointer; }
  constexpr index size() const noexcept { return N; }

  // Element Access

  constexpr reference access(index n) const noexcept {
    assert(0 <= n && n < N);
    return *(m_pointer + n);
  }

  constexpr reference operator[](index n) const noexcept { return access(n); }

  constexpr reference operator()(index n) const noexcept { return access(n); }
  constexpr reference operator()(const array<index, 1>& n) const noexcept {
    return access(n[0]);
  }

  // Iterators

  constexpr iterator begin() const noexcept { return m_pointer; }
  constexpr const_iterator cbegin() const noexcept { return m_pointer; }

  constexpr iterator end() const noexcept { return m_pointer + N; }
  constexpr const_iterator cend() const noexcept { return m_pointer + N; }

  // Returns true if this span points to something different than nullptr.
  constexpr operator bool() const noexcept { return m_pointer != nullptr; }

private:
  pointer m_pointer{nullptr};
};

/// Class specialisation for the dynamically sized case.
template <typename T> class span<T, dyn> {
public:
  // Types

  using element_type = T;
  using pointer = T*;
  using reference = T&;
  using iterator = pointer;
  using const_iterator = const T*;

  // Constructors

  constexpr span() = default;

  template <
      typename S, index M,
      typename = std::enable_if_t<std::is_convertible<S*, pointer>::value>>
  constexpr span(const span<S, M>& other) noexcept
      : m_pointer{other.data()}, m_size{other.size()} {}

  template <
      typename S, index M,
      typename = std::enable_if_t<std::is_convertible<S*, pointer>::value>>
  constexpr span operator=(const span<S, M>& other) noexcept {
    m_pointer = other.data();
    m_size = other.size();
    return *this;
  }

  template <std::size_t N>
  constexpr span(element_type (&a)[N]) noexcept
      : m_pointer{a.data()}, m_size{static_cast<index>(N)} {}

  constexpr span(pointer ptr, index size) noexcept
      : m_pointer{ptr}, m_size{size} {}

  template <std::size_t N, typename Alloc, typename S = element_type,
            typename = std::enable_if_t<std::is_const<S>::value>>
  constexpr span(const array<std::remove_const_t<element_type>, N>& a) noexcept
      : m_pointer{a.data()}, m_size{static_cast<index>(N)} {}

  template <std::size_t N>
  constexpr span(array<std::remove_const_t<element_type>, N>& a) noexcept
      : m_pointer{a.data()}, m_size{static_cast<index>(N)} {}

  template <typename Alloc>
  constexpr span(
      std::vector<std::remove_const_t<element_type>, Alloc>& v) noexcept
      : m_pointer{v.data()}, m_size{static_cast<index>(v.size())} {}

  template <typename Alloc, typename S = element_type,
            typename = std::enable_if_t<std::is_const<S>::value>>
  constexpr span(
      const std::vector<std::remove_const_t<element_type>, Alloc>& v) noexcept
      : m_pointer{v.data()}, m_size{static_cast<index>(v.size())} {}

  // Class Member Access

  constexpr pointer data() const noexcept { return m_pointer; }
  constexpr index size() const noexcept { return m_size; }

  // Element Access

  constexpr reference access(index n) const noexcept {
    assert(0 <= n && n < m_size);
    return *(m_pointer + n);
  }

  constexpr reference operator[](index n) const noexcept { return access(n); }

  constexpr reference operator()(index n) const noexcept { return access(n); }
  constexpr reference operator()(const array<index, 1>& n) const noexcept {
    return access(n[0]);
  }

  // Iterators

  constexpr iterator begin() const noexcept { return m_pointer; }
  constexpr const_iterator cbegin() const noexcept { return m_pointer; }

  constexpr iterator end() const noexcept { return m_pointer + m_size; }
  constexpr const_iterator cend() const noexcept { return m_pointer + m_size; }

  // Returns true if this span points to something different than nullptr.
  constexpr operator bool() const noexcept { return m_pointer != nullptr; }

private:
  pointer m_pointer{nullptr};
  index m_size{0};
};

template <typename T> struct is_span : std::false_type {};
template <typename T, index N> struct is_span<span<T, N>> : std::true_type {};
template <typename T> static constexpr bool is_span_v = is_span<T>::value;

////////////////////////////////////////////////////////////////////////////////
// make_span

template <typename T, std::size_t N>
auto make_span(std::array<T, N>& array) noexcept {
  return span<T, N>(array);
}

template <typename T, std::size_t N>
auto make_span(const std::array<T, N>& array) noexcept {
  return span<const T, N>(array);
}

template <typename T, typename Alloc>
auto make_span(std::vector<T, Alloc>& array) noexcept {
  return span<T>(array);
}

template <typename T, typename Alloc>
auto make_span(const std::vector<T, Alloc>& array) noexcept {
  return span<const T>(array);
}

////////////////////////////////////////////////////////////////////////////////
// drop / take with static extents

template <index N, typename T, index Size>
constexpr std::enable_if_t<(0 <= N && N < Size), span<T, Size - N>>
drop(span<T, Size> view) noexcept {
  return span<T, Size - N>(view.data() + N, Size - N);
}

template <index N, typename T, index Size>
constexpr std::enable_if_t<(0 < N && N <= Size), span<T, N>>
take(span<T, Size> view) noexcept {
  return span<T, N>(view.data(), N);
}

////////////////////////////////////////////////////////////////////////////////
// drop / take with dynamic extents

template <index N, typename T> span<T> drop(span<T> view) {
  return span<T>(view.data() + std::min(view.size(), N),
                 std::max(0, view.size() - N));
}

template <typename T>
span<T, dyn> drop(const span<T, dyn>& view, index n) noexcept {
  return span<T>(view.data() + std::min(view.size(), n),
                 std::max(0, view.size() - n));
}

template <index N, typename T> span<T, N> take(span<T> view) {
  if (view.size() < N) {
    throw std::out_of_range{"take: span's length is smaller than than."};
  }
  return span<T, N>(view.data());
}

} // namespace fub

#endif // !SPAN_HPP
