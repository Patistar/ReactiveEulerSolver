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

#ifndef FUB_GRID_VARIABLE_REF_HPP
#define FUB_GRID_VARIABLE_REF_HPP

// #include "fub/accessor_simd.hpp"
#include "fub/variable_list.hpp"

namespace fub {
inline namespace v1 {
template <typename VariableList, typename T,
          typename Accessor = accessor_native<T>>
class variable_ref : VariableList, Accessor {
public:
  using accessor = Accessor;
  using value_type = typename Accessor::value_type;
  using pointer = typename Accessor::pointer;
  using reference = typename Accessor::reference;

  variable_ref() = delete;
  variable_ref(const variable_ref&) = delete;
  variable_ref(variable_ref&&) = default;

  template <typename S, typename OtherAccessor,
            typename = std::enable_if_t<!std::is_same<T, S>{}>,
            typename = std::enable_if_t<std::is_convertible<
                typename OtherAccessor::pointer, pointer>{}>>
  variable_ref(variable_ref<VariableList, S, OtherAccessor>&& other)
      : VariableList(other.get_variable_list()), m_pointer{other.get_pointer()},
        m_stride{other.get_stride()} {}

  variable_ref(VariableList list, pointer ptr, std::ptrdiff_t stride,
               Accessor accessor = Accessor()) noexcept
      : VariableList(list), Accessor(accessor), m_pointer(ptr),
        m_stride(stride) {}

  template <typename Variables>
  variable_ref& operator=(const Variables& variables) {
    for_each(get_variable_list(), [&](auto variable) {
      this->operator[](variable) = variables[variable];
    });
    return *this;
  }

  template <typename Tag>
  constexpr reference operator[](Tag tag) const noexcept {
    std::ptrdiff_t size = get_variable_list().size() * m_stride;
    std::ptrdiff_t index = *get_variable_list().index(tag) * m_stride;
    return get_accessor().access(m_pointer, size, index);
  }

  constexpr VariableList get_variable_list() const noexcept { return *this; }
  constexpr Accessor get_accessor() const noexcept { return *this; }
  constexpr std::ptrdiff_t get_stride() const noexcept { return m_stride; }
  constexpr pointer get_pointer() const noexcept { return m_pointer; }

private:
  pointer m_pointer;
  std::ptrdiff_t m_stride;
};

template <typename VariableList, typename T,
          typename Accessor = accessor_native<T>>
class variable_iterator : VariableList, Accessor {
public:
  using pointer = typename Accessor::pointer;
  using reference = variable_ref<VariableList, T, Accessor>;
  using value_type = reference;
  using difference_type = std::ptrdiff_t;
  using iterator_category = std::forward_iterator_tag;

  variable_iterator() = default;

  variable_iterator(VariableList list, pointer ptr, std::ptrdiff_t stride,
                    Accessor accessor = Accessor()) noexcept
      : VariableList(list), Accessor(accessor), m_pointer(ptr),
        m_stride(stride) {}

  variable_ref<VariableList, T, Accessor> read() const noexcept {
    return {get_variable_list(), m_pointer, m_stride, get_accessor()};
  }

  constexpr VariableList get_variable_list() const noexcept { return *this; }
  constexpr Accessor get_accessor() const noexcept { return *this; }
  constexpr pointer get_pointer() const noexcept { return m_pointer; }
  constexpr std::ptrdiff_t get_stride() const noexcept { return m_stride; }

  variable_iterator& operator++() noexcept {
    ++m_pointer;
    return *this;
  }

  variable_iterator& operator--() noexcept {
    --m_pointer;
    return *this;
  }

  variable_iterator operator++(int) noexcept {
    variable_iterator old{*this};
    ++m_pointer;
    return old;
  }

  variable_iterator& operator--(int) noexcept {
    variable_iterator old{*this};
    --m_pointer;
    return old;
  }

  bool operator==(const variable_iterator& other) const noexcept {
    return get_variable_list() == other.get_variable_list() &&
           get_accessor() == other.get_accessor() &&
           m_pointer == other.m_pointer;
  }

  bool operator!=(const variable_iterator& other) const noexcept {
    return !(*this == other);
  }

  reference operator*() const noexcept {
    return reference{get_variable_list(), m_pointer, m_stride, get_accessor()};
  }

private:
  pointer m_pointer{nullptr};
  std::ptrdiff_t m_stride{};
};

// template <typename VariableList, typename T, typename Accessor>
// variable_iterator<
//     VariableList, T,
//     accessor_simd<T, simd_abi::native<T>, flags::vector_aligned_tag>>
// simdify(variable_iterator<VariableList, T, Accessor> it) {
//   accessor_simd<T, simd_abi::native<T>, flags::vector_aligned_tag> a{};
//   return {it.get_variable_list(), it.get_pointer(),
//           accessible_size(a, it.get_pointer(), it.get_stride())};
// }

} // namespace v1
} // namespace fub

#endif