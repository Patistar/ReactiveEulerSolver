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
#include "fub/native_accessor.hpp"
#include "fub/type_traits.hpp"

namespace fub {
template <typename T, index N, typename Accessor>
class basic_span_storage : private Accessor {
public:
  static_assert(N > 0, "N must be positive for statically sized spans.");

  using pointer = typename Accessor::template pointer<T>;
  using const_pointer = typename Accessor::template pointer<const T>;
  using reference = typename Accessor::template reference<T>;

  constexpr basic_span_storage() = delete;

  constexpr explicit basic_span_storage(pointer p,
                                        const Accessor& a = Accessor())
      : Accessor(a), m_pointer{p} {}

  index size() const noexcept { return N; }

  reference access(index n) const
      noexcept(noexcept(Accessor::access(m_pointer, n))) {
    return Accessor::access(m_pointer, n);
  }

  pointer get_pointer() const noexcept { return m_pointer; }

private:
  pointer m_pointer;
};

template <typename T, typename Accessor>
class basic_span_storage<T, dyn, Accessor> : private Accessor {
public:
  using pointer = typename Accessor::template pointer<T>;
  using const_pointer = typename Accessor::template pointer<const T>;
  using reference = typename Accessor::template reference<T>;

  constexpr basic_span_storage() = default;

  constexpr explicit basic_span_storage(pointer p, index sz,
                                        const Accessor& a = Accessor())
      : Accessor(a), m_pointer{p}, m_size{sz} {}

  index size() const noexcept { return m_size; }

  reference access(index n) const { return Accessor::access(m_pointer, n); }

  pointer get_pointer() const noexcept { return m_pointer; }

private:
  pointer m_pointer{nullptr};
  index m_size{0};
};

template <typename R> using data_t = decltype(fub::data(std::declval<R>()));

template <typename A> struct is_std_array : bool_constant<false> {};

template <typename T, std::size_t N>
struct is_std_array<array<T, N>> : bool_constant<true> {};

template <typename T, index N, typename Accessor> class basic_span;

template <typename T> struct is_span : std::false_type {};

template <typename T, index N, typename A>
struct is_span<basic_span<T, N, A>> : std::true_type {};

template <typename T> static constexpr bool is_span_v = is_span<T>::value;

template <typename T, index N, typename Accessor> class basic_span {
public:
  using element_type = T;
  using storage_type = basic_span_storage<T, N, Accessor>;
  using value_type = std::remove_cv_t<T>;
  using pointer = typename storage_type::pointer;
  using const_pointer = typename storage_type::const_pointer;
  using reference = typename storage_type::reference;
  using iterator = pointer;
  using const_iterator = const_pointer;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;
  using index_type = index;

  static constexpr index_type extent = N;

  /////////////////////////////////////////////////////////////////////////////
  // Constructors

  constexpr basic_span(pointer p, [[maybe_unused]] index_type size) noexcept
      : m_storage{p} {
    assert(N <= size);
  }
  constexpr basic_span(pointer first, [[maybe_unused]] pointer last) noexcept
      : m_storage{first} {
    assert(N <= (last - first));
  }

  template <std::size_t M, typename = std::enable_if_t<N <= M>>
  constexpr basic_span(element_type (&arr)[M]) noexcept // NOLINT
      : m_storage{&arr[0]} {}

  template <typename S, std::size_t M,
            typename = std::enable_if_t<std::is_convertible<
                typename array<S, M>::const_pointer, pointer>::value>,
            typename = std::enable_if_t<N <= M>>
  constexpr basic_span(const array<S, M>& arr) noexcept // NOLINT
      : m_storage{arr.data()} {}

  template <typename S, std::size_t M,
            typename = std::enable_if_t<std::is_convertible<
                typename array<S, M>::pointer, pointer>::value>,
            typename = std::enable_if_t<N <= M>>
  constexpr basic_span(array<S, M>& arr) noexcept // NOLINT
      : m_storage{arr.data()} {}

  template <typename S, index_type M, typename A,
            typename = std::enable_if_t<std::is_convertible<
                typename basic_span<S, M, A>::pointer, pointer>::value>,
            typename = std::enable_if_t<N <= M>>
  constexpr basic_span(const basic_span<S, M, A>& s) noexcept // NOLINT
      : m_storage{s.data()} {}

  constexpr basic_span(const basic_span& s) noexcept = default;

  template <
      typename Container,
      typename =
          std::enable_if_t<!is_std_array<std::decay_t<Container>>::value>,
      typename = std::enable_if_t<!is_span<std::decay_t<Container>>::value>,
      typename = std::enable_if_t<
          std::is_convertible<detected_t<data_t, Container>, pointer>::value>>
  constexpr basic_span(Container&& container) noexcept // NOLINT
      : m_storage{fub::data(container)} {
    assert(N <= fub::size(container));
  }

  /////////////////////////////////////////////////////////////////////////////
  // Element Access

  constexpr reference operator[](index_type n) const {
    return m_storage.access(n);
  }

  constexpr reference operator()(index_type n) const {
    return m_storage.access(n);
  }

  /////////////////////////////////////////////////////////////////////////////
  // Member Variable Accessors

  constexpr index_type size() const noexcept { return m_storage.size(); }

  constexpr pointer data() const noexcept { return m_storage.get_pointer(); }

  /////////////////////////////////////////////////////////////////////////////
  // Iterators

  constexpr iterator begin() const noexcept { return m_storage.get_pointer(); }
  constexpr const_iterator cbegin() const noexcept {
    return m_storage.get_pointer();
  }
  constexpr reverse_iterator rbegin() const noexcept {
    return {m_storage.get_pointer()};
  }
  constexpr const_reverse_iterator crbegin() const noexcept {
    return {m_storage.get_pointer()};
  }

  constexpr iterator end() const noexcept {
    return m_storage.get_pointer() + size();
  }
  constexpr iterator cend() const noexcept {
    return m_storage.get_pointer() + size();
  }
  constexpr iterator rend() const noexcept {
    return {m_storage.get_pointer() + size()};
  }
  constexpr iterator crend() const noexcept {
    return {m_storage.get_pointer() + size()};
  }

  /// Returns always true since this span version can not be emtpy.
  ///
  /// @note This makes span<T, N> contextually convertible to bool.
  ///
  /// @example
  ///     span<int, 3> s = ...;
  ///     // ...
  ///     if (s) {
  ///         // s points to some valid array
  ///     } else {
  ///         // s is empty
  ///     }
  constexpr explicit operator bool() const noexcept { return true; }

private:
  storage_type m_storage;
};

template <typename T, typename Accessor> class basic_span<T, dyn, Accessor> {
public:
  using element_type = T;
  using storage_type = basic_span_storage<T, dyn, Accessor>;
  using value_type = std::remove_cv_t<T>;
  using pointer = typename storage_type::pointer;
  using const_pointer = typename storage_type::const_pointer;
  using reference = typename storage_type::reference;
  using iterator = pointer;
  using const_iterator = const_pointer;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;
  using index_type = index;

  static constexpr index_type extent = dyn;

  /////////////////////////////////////////////////////////////////////////////
  // Constructors

  constexpr basic_span() = default;

  constexpr basic_span(pointer p, index_type size) noexcept
      : m_storage{p, size} {}
  constexpr basic_span(pointer first, pointer last) noexcept
      : m_storage{first, last - first} {}

  template <std::size_t M>
  constexpr basic_span(element_type (&arr)[M]) noexcept // NOLINT
      : m_storage{&arr[0], static_cast<index_type>(M)} {}

  template <typename S, std::size_t M,
            typename = std::enable_if_t<std::is_convertible<
                typename array<S, M>::const_pointer, pointer>::value>>
  constexpr basic_span(const array<S, M>& arr) noexcept // NOLINT
      : m_storage{arr.data(), arr.size()} {}

  template <typename S, std::size_t M,
            typename = std::enable_if_t<std::is_convertible<
                typename array<S, M>::pointer, pointer>::value>>
  constexpr basic_span(array<S, M>& arr) noexcept // NOLINT
      : m_storage{arr.data(), arr.size()} {}

  template <
      typename Container,
      typename =
          std::enable_if_t<!is_std_array<std::decay_t<Container>>::value>,
      typename = std::enable_if_t<!is_span<std::decay_t<Container>>::value>,
      typename = std::enable_if_t<
          std::is_convertible<detected_t<data_t, Container>, pointer>::value>>
  constexpr basic_span(Container&& container) noexcept // NOLINT
      : m_storage{container.data(),
                  static_cast<std::ptrdiff_t>(container.size())} {}

  template <typename S, index_type M, typename A,
            typename = std::enable_if_t<std::is_convertible<
                typename basic_span<S, M, A>::pointer, pointer>::value>>
  constexpr basic_span(const basic_span<S, M, A>& s) noexcept // NOLINT
      : m_storage{s.data(), s.size()} {}

  constexpr basic_span(const basic_span& s) noexcept = default;

  /////////////////////////////////////////////////////////////////////////////
  // Element Access

  constexpr reference operator[](index_type n) const {
    return m_storage.access(n);
  }

  constexpr reference operator()(index_type n) const {
    return m_storage.access(n);
  }

  /////////////////////////////////////////////////////////////////////////////
  // Member Variable Accessors

  constexpr index_type size() const noexcept { return m_storage.size(); }

  constexpr pointer data() const noexcept { return m_storage.get_pointer(); }

  /////////////////////////////////////////////////////////////////////////////
  // Iterators

  constexpr iterator begin() const noexcept { return m_storage.get_pointer(); }
  constexpr const_iterator cbegin() const noexcept {
    return m_storage.get_pointer();
  }
  constexpr reverse_iterator rbegin() const noexcept {
    return {m_storage.get_pointer()};
  }
  constexpr const_reverse_iterator crbegin() const noexcept {
    return {m_storage.get_pointer()};
  }

  constexpr iterator end() const noexcept {
    return m_storage.get_pointer() + size();
  }
  constexpr iterator cend() const noexcept {
    return m_storage.get_pointer() + size();
  }
  constexpr iterator rend() const noexcept {
    return {m_storage.get_pointer() + size()};
  }
  constexpr iterator crend() const noexcept {
    return {m_storage.get_pointer() + size()};
  }

  /// Returns true if this span points to something different than nullptr.
  ///
  /// @note This makes span<T> contextually convertible to bool.
  ///
  /// @example
  ///     span<int> s = nullptr;
  ///     // ...
  ///     if (s) {
  ///         // s points to some valid array
  ///     } else {
  ///         // s is empty
  ///     }
  constexpr explicit operator bool() const noexcept {
    return data() != nullptr;
  }

private:
  storage_type m_storage;
};

template <typename T, index Extent = dyn>
using span = basic_span<T, Extent, native_accessor>;

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
                 std::max(index{}, view.size() - N));
}

template <typename T> span<T> drop(const span<T>& view, index n) noexcept {
  return span<T>(view.data() + std::min(view.size(), n),
                 std::max(index{}, view.size() - n));
}

template <index N, typename T> span<T, N> take(span<T> view) {
  if (view.size() < N) {
    throw std::out_of_range{"take: span's length is smaller than than."};
  }
  return span<T, N>(view.data(), N);
}

} // namespace fub

#endif // !SPAN_HPP
