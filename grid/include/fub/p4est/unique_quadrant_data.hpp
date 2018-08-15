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

#ifndef FUB_P4EST_UNIQUE_QUADRANT_DATA_HPP
#define FUB_P4EST_UNIQUE_QUADRANT_DATA_HPP

#include "fub/p4est/quadrant.hpp"
#include "fub/span.hpp"

#include <Vc/vector.h>

namespace fub {
inline namespace v1 {
namespace p4est {

template <typename T> class unique_quadrant_data;

template <typename T> class unique_quadrant_data<T[]> {
public:
  static_assert(std::is_trivially_destructible<T>::value,
                "T is not trivially destructible and can not be used with "
                "unique_quadrant_data.");
  static_assert(std::is_trivially_copy_assignable<T>::value,
                "T is not trivially copy assignable and can not be used with "
                "unique_quadrant_data.");

  using allocator = boost::container::pmr::polymorphic_allocator<T>;

  template <typename S>
  using rebind_alloc =
      typename std::allocator_traits<allocator>::template rebind_alloc<S>;

  /// Constructs an empty quadrant data object.
  unique_quadrant_data() = default;

  /// Constructs data arrays of size `size` for each quadrant in `quadrants`.
  unique_quadrant_data(std::ptrdiff_t num_quadrants, std::ptrdiff_t size,
                       allocator alloc = allocator())
      : m_allocator(alloc), m_quadrant_to_data{} {
    if (num_quadrants > 0) {
      m_quadrant_to_data = span<span<T>>(
          rebind_alloc<span<T>>(m_allocator).allocate(num_quadrants),
          num_quadrants);
    }
    for (fub::span<T>& span : m_quadrant_to_data) {
      new (&span) fub::span<T>{m_allocator.allocate(size), size};
    }
  }

  unique_quadrant_data(const unique_quadrant_data&) = delete;
  unique_quadrant_data& operator=(const unique_quadrant_data&) = delete;

  /// Move ownership of the data from other to this.
  unique_quadrant_data(unique_quadrant_data&& other) noexcept
      : m_allocator{other.get_allocator()}, m_quadrant_to_data{
                                                other.release()} {}

  /// Move ownership of the data from other to this.
  unique_quadrant_data& operator=(unique_quadrant_data&& other) noexcept {
    reset(other.release(), other.get_allocator());
    return *this;
  }

  ~unique_quadrant_data() noexcept {
    if (m_quadrant_to_data) {
      for (fub::span<T> span : m_quadrant_to_data) {
        m_allocator.deallocate(span.data(), span.size());
      }
      rebind_alloc<span<T>>(m_allocator)
          .deallocate(m_quadrant_to_data.data(), m_quadrant_to_data.size());
    }
  }

  std::ptrdiff_t size() const noexcept { return m_quadrant_to_data.size(); }

  span<const T> operator[](std::ptrdiff_t idx) const noexcept {
    assert(0 <= idx && idx < size());
    return m_quadrant_to_data[idx];
  }

  span<T> operator[](std::ptrdiff_t idx) noexcept {
    assert(0 <= idx && idx < size());
    return m_quadrant_to_data[idx];
  }

  allocator get_allocator() const noexcept { return m_allocator; }

  span<span<T>> release() noexcept {
    return std::exchange(m_quadrant_to_data, span<span<T>>());
  }

  void reset(span<span<T>> quad_to_data, allocator alloc) noexcept {
    this->~unique_quadrant_data();
    m_allocator = alloc;
    m_quadrant_to_data = quad_to_data;
  }

  span<const span<T>> data() noexcept { return m_quadrant_to_data; }

  span<const span<const T>> data() const noexcept {
    return span_cast<const span<const T>>(m_quadrant_to_data);
  }

private:
  allocator m_allocator{};
  span<span<T>> m_quadrant_to_data{};
};

} // namespace p4est
} // namespace v1
} // namespace fub

#endif