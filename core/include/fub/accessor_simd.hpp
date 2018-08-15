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

#ifndef FUB_CORE_ACCESSOR_SIMD_HPP
#define FUB_CORE_ACCESSOR_SIMD_HPP

#include "fub/aligned_address.hpp"
#include "fub/fancy_pointer.hpp"

#include <Vc/vector.h>

#include <cassert>

namespace fub {
inline namespace v1 {
template <typename T, typename Abi, typename AlignmentFlag>
struct simd_pointer_helper {
  using type = T*;
};

template <typename T, typename Abi>
struct simd_pointer_helper<T, Abi, Vc::AlignedTag> {
  using simd_type = Vc::Vector<remove_cvref_t<T>, Abi>;
  using type = fancy_pointer<T, aligned_address<simd_type::MemoryAlignment>>;
};

template <typename T, typename Abi, typename Alignment>
using simd_pointer = typename simd_pointer_helper<T, Abi, Alignment>::type;

template <typename T, typename Abi, typename Alignment> class simd_reference {
public:
  using element_type = T;
  using value_type = remove_cvref_t<T>;
  using pointer = simd_pointer<T, Abi, Alignment>;
  using simd_type = Vc::Vector<value_type, Abi>;

  simd_reference(pointer p) noexcept : m_pointer{p} {}
  simd_reference(const simd_reference&) = delete;
  simd_reference& operator=(const simd_reference&) = delete;
  simd_reference(simd_reference&&) = default;
  simd_reference& operator=(simd_reference&&) = default;

  operator simd_type() const noexcept {
    return simd_type{static_cast<T*>(m_pointer), Alignment()};
  }

  simd_reference& operator=(const simd_type& v) noexcept {
    v.store(static_cast<T*>(m_pointer), Alignment());
    return *this;
  }

private:
  pointer m_pointer;
};

template <typename T, typename Abi, typename Alignment>
class simd_reference<const T, Abi, Alignment> {
public:
  using value_type = T;
  using element_type = const T;
  using pointer = simd_pointer<const T, Abi, Alignment>;
  using simd_type = Vc::Vector<remove_cvref_t<T>, Abi>;

  simd_reference(pointer p) noexcept : m_pointer{p} {}
  simd_reference(const simd_reference&) = delete;
  simd_reference& operator=(const simd_reference&) = delete;
  simd_reference(simd_reference&&) = default;
  simd_reference& operator=(simd_reference&&) = default;

  operator simd_type() const noexcept {
    return simd_type{static_cast<const T*>(m_pointer), Alignment()};
  }

private:
  pointer m_pointer;
};

template <typename T, typename Abi = Vc::VectorAbi::Best<T>, typename Alignment = Vc::AlignedTag> struct accessor_simd {
  using element_type = T;
  using simd_type = Vc::Vector<remove_cvref_t<T>, Abi>;
  using value_type = simd_type;
  using pointer = simd_pointer<T, Abi, Alignment>;
  using reference = simd_reference<T, Abi, Alignment>;

  static constexpr int alignment() noexcept {
    return simd_type::MemoryAlignment;
  }

  static constexpr pointer to_pointer(element_type& element) {
    return static_cast<void*>(&element);
  }

  static constexpr reference access(pointer origin, std::ptrdiff_t /* size */,
                                    std::ptrdiff_t offset) noexcept {
    pointer ptr = origin + offset;
    return reference{ptr};
  }

  template <typename S> using rebind = accessor_simd<S, Abi, Vc::AlignedTag>;
};

template <typename T, typename Abi, typename Alignment>
constexpr bool operator==(accessor_simd<T, Abi, Alignment>,
                          accessor_simd<T, Abi, Alignment>) noexcept {
  return true;
}

template <typename T, typename Abi, typename Alignment>
constexpr bool operator!=(accessor_simd<T, Abi, Alignment>,
                          accessor_simd<T, Abi, Alignment>) noexcept {
  return false;
}

} // namespace v1
} // namespace fub

#endif
