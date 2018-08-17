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

#include <p4est_communication.h>

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

  /// Constructs an empty quadrant data object.
  unique_quadrant_data() = default;

  /// Constructs data arrays of size `size` for each quadrant in `quadrants`.
  unique_quadrant_data(std::ptrdiff_t num_quadrants, std::ptrdiff_t size,
                       allocator alloc = allocator())
      : m_stride{size}, m_data(num_quadrants * size,
                               std::numeric_limits<T>::signaling_NaN(), alloc) {
  }

  std::ptrdiff_t data_size() const noexcept { return m_stride; }

  std::ptrdiff_t size() const noexcept { return m_data.size() / m_stride; }

  fub::span<const T> operator[](std::ptrdiff_t idx) const noexcept {
    assert(0 <= idx && idx < size());
    return {&m_data[idx * m_stride], m_stride};
  }

  fub::span<T> operator[](std::ptrdiff_t idx) noexcept {
    assert(0 <= idx && idx < size());
    return {&m_data[idx * m_stride], m_stride};
  }

  allocator get_allocator() const noexcept { return m_data.get_allocator(); }

  fub::span<T> span() noexcept { return m_data; }

  fub::span<const T> span() const noexcept { return m_data; }

private:
  std::ptrdiff_t m_stride{};
  std::vector<T, allocator> m_data{};
};

template <typename T>
unique_quadrant_data<T[]>
transfer_data(const forest<2>& new_forest, const forest<2>& old_forest,
              const unique_quadrant_data<T[]>& old_data) {
  span<const p4est_gloidx_t> dest_gfp = new_forest.global_first_position();
  span<const p4est_gloidx_t> src_gfp = old_forest.global_first_position();
  MPI_Comm comm = new_forest.mpi_communicator();
  std::ptrdiff_t data_size = old_data.data_size();
  unique_quadrant_data<T[]> new_data(new_forest.local_num_quadrants(),
                                     data_size);
  p4est_transfer_fixed(dest_gfp.data(), src_gfp.data(), comm, 0, new_data.span().data(),
                       old_data.span().data(), sizeof(T) * data_size);
  return new_data;
}

} // namespace p4est
} // namespace v1
} // namespace fub

#endif