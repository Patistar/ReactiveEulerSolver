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

#ifndef FUB_CORE_LAYOUT_RIGHT_HPP
#define FUB_CORE_LAYOUT_RIGHT_HPP

#include "fub/algorithm.hpp"
#include "fub/optional.hpp"
#include "fub/type_traits.hpp"

#include <array>
#include <functional>

namespace fub {

struct layout_right {
  template <typename Extents> class mapping : private Extents {
  public:
    using index_type = typename Extents::index_type;

    constexpr mapping() = default;
    constexpr mapping(const mapping&) = default;
    constexpr mapping(mapping&&) = default;
    mapping& operator=(const mapping&) noexcept = default;
    mapping& operator=(mapping&&) noexcept = default;

    /// \note This is explicit opposed to P0009
    explicit constexpr mapping(const Extents& extents) : Extents{extents} {}

    constexpr const Extents& get_extents() const noexcept { return *this; }

    constexpr index_type required_span_size() const noexcept {
      return size(get_extents());
    }

    template <
        typename... IndexType,
        typename = std::enable_if_t<
            conjunction<std::is_convertible<IndexType, index_type>...>::value>,
        typename = std::enable_if_t<(sizeof...(IndexType) == Extents::rank())>>
    constexpr index_type operator()(IndexType... indices) const noexcept {
      return do_mapping(std::make_index_sequence<Extents::rank()>(),
                        indices...);
    }

    constexpr index_type stride(std::size_t r) const noexcept {
      index_type stride = 1;
      for (std::size_t dim = r + 1; dim < Extents::rank(); ++dim) {
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

template <typename Extents>
constexpr typename Extents::index_type
static_required_span_size(layout_right, const Extents&) noexcept {
  typename Extents::index_type sz = size(Extents());
  return sz ? sz : dynamic_extent;
}


template <typename Extents>
constexpr index_array_t<Extents>
next(const layout_right::mapping<Extents>& mapping,
     index_array_t<Extents> index) noexcept {
  assert(is_in_range(mapping.get_extents(), index));
  std::size_t r = Extents::rank();
  index[r - 1] += 1;
  while (0 < r && index[r - 1] >= mapping.get_extents().extent(r - 1)) {
    index[r - 1] = 0;
    r -= 1;
    index[r] += 1;
  }
  assert(is_in_range(mapping.get_extents(), index));
  return index;
}

template <typename Extents, typename Function>
constexpr Function for_each_index(const layout_right::mapping<Extents>& mapping,
                                  Function fun) {
  using index_type = typename Extents::index_type;
  Extents extents = mapping.get_extents();
  index_type sz = size(extents);
  if (sz == 0) {
    return fun;
  }
  std::array<index_type, Extents::rank()> index{};
  while (sz > 0) {
    fub::invoke(fun, index);
    index = next(mapping, index);
    sz -= 1;
  }
  return fun;
}

} // namespace fub

#endif
