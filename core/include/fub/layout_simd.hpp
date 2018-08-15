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

#ifndef FUB_CORE_LAYOUT_SIMD_HPP
#define FUB_CORE_LAYOUT_SIMD_HPP

#include "fub/algorithm.hpp"
#include "fub/extents.hpp"
#include "fub/simd.hpp"
#include "fub/type_traits.hpp"

#include <Vc/vector.h>

#include <array>

namespace fub {
inline namespace v1 {

template <typename T, typename Abi> struct layout_simd_padded {
  template <typename Extents> class mapping : private Extents {
  public:
    static_assert(is_extents_v<Extents>,
                  "Extents must satisfy the extents concept.");

    using index_type = typename Extents::index_type;
    using value_type = remove_cvref_t<T>;
    using simd_type = Vc::Vector<value_type, Abi>;

    // CONSTRUCTORS

    constexpr mapping() = default;
    constexpr mapping(const mapping&) = default;
    constexpr mapping(mapping&&) = default;
    mapping& operator=(const mapping& other) noexcept = default;
    mapping& operator=(mapping&& other) noexcept = default;

    /// Implicit Conversion from Extents.
    ///
    /// Throws: Nothing.
    ///
    /// Postcondition: get_extents() == extents.
    constexpr mapping(const Extents& extents) : Extents(extents) {}

    // OBSERVERS

    constexpr const Extents& get_extents() const noexcept { return *this; }

    /// Returns the required span size after adding padding in each dimension
    /// according to the specified simd alignment.
    ///
    /// Throws: Nothing.
    constexpr index_type required_span_size() const noexcept {
      constexpr index_type width = simd_type::MemoryAlignment;
      index_type factor = 1;
      for (std::size_t dim = 0; dim < Extents::rank(); ++dim) {
        factor *= width;
      }
      return size(get_extents()) * factor;
    }

    /// Returns the codomain index for specified multi dimensional index
    /// coordinates.
    ///
    /// Effect: Equivalent to offset in
    ///
    ///    index_type offset = 0;
    ///    for(size_t k = 0; k < Extents::rank(); ++k) {
    ///      offset += i[k] * stride(k);
    ///    }
    ///
    /// Throws: Nothing.
    template <
        typename... IndexType,
        typename = std::enable_if_t<
            conjunction<std::is_convertible<IndexType, index_type>...>::value>,
        typename = std::enable_if_t<(sizeof...(IndexType) == Extents::rank())>>
    constexpr index_type operator()(IndexType... indices) const noexcept {
      return do_mapping(std::make_index_sequence<Extents::rank()>(),
                        indices...);
    }

    /// Returns the stride value in a specified dimension r.
    ///
    /// Throws: Nothing.
    constexpr index_type stride(std::size_t r) const noexcept {
      index_type stride = 1;
      for (std::size_t dim = 0; dim < r; ++dim) {
        stride *= get_extents().extent(dim);
      }
      return stride;
    }

    static constexpr bool is_always_unique() noexcept { return true; }
    static constexpr bool is_always_contiguous() noexcept { return true; }
    static constexpr bool is_always_strided() noexcept { return true; }

    constexpr bool is_unique() const noexcept { return true; }
    constexpr bool is_contiguous() const noexcept { return true; }
    constexpr bool is_strided() const noexcept { return true; }

    template <class OtherExtents>
    constexpr bool operator==(const mapping<OtherExtents>& other) const
        noexcept {
      return get_extents() == other.get_extents();
    }

    template <class OtherExtents>
    constexpr bool operator!=(const mapping<OtherExtents>& other) const
        noexcept {
      return !(*this == other);
    }

  private:
    template <std::size_t... Is, typename... IndexType>
    constexpr index_type do_mapping(std::index_sequence<Is...>,
                                    IndexType... indices) const noexcept {
      const index_type is[Extents::rank()]{indices...};
      const index_type strides[Extents::rank()]{stride(Is)...};
      index_type sum = fub::transform_reduce(strides, is, index_type(0));
      return sum;
    }
  };
};

template <typename T, typename Abi, typename Extents>
constexpr typename Extents::index_type
static_required_span_size(layout_simd_padded<T, Abi>,
                          const Extents& extents) noexcept {
  typename Extents::index_type sz = size(Extents());
  if (sz) {
    typename layout_simd_padded<T, Abi>::template mapping<Extents> m(extents);
    return m.required_span_size();
  }
  return dynamic_extent;
}

} // namespace v1
} // namespace fub

#endif
