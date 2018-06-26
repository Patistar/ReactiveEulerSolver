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

#ifndef FUB_CORE_LAYOUT_LEFT_HPP
#define FUB_CORE_LAYOUT_LEFT_HPP

#include "fub/algorithm.hpp"
#include "fub/optional.hpp"
#include "fub/type_traits.hpp"

#include <array>
#include <functional>

namespace fub {

struct layout_left {
  template <typename Extents> class mapping : private Extents {
  public:
    static constexpr int rank = Extents::rank;

    using Extents::Extents;

    explicit constexpr mapping(const Extents& extents) : Extents{extents} {}

    constexpr const Extents& extents() const noexcept { return *this; }

    template <typename... Indices,
              typename = std::enable_if_t<(sizeof...(Indices) == rank)>>
    constexpr index operator()(Indices... indices) const noexcept {
      return do_mapping(std::array<index, rank>{{indices...}},
                        make_int_sequence<rank>());
    }

    template <int R> constexpr index stride() const noexcept {
      return do_stride(make_int_sequence<R>());
    }

    constexpr optional<std::array<index, rank>>
    next(const std::array<index, rank>& indices) const noexcept {
      return do_next(indices, int_c<0>);
    }

  private:
    template <int Dim>
    constexpr optional<std::array<index, rank>>
    do_next(std::array<index, rank> indices, int_constant<Dim>) const noexcept {
      static_assert(0 <= Dim && Dim < Extents::rank,
                    "Dimension is out of range.");
      ++indices[Dim];
      if (indices[Dim] < extents().get(Dim)) {
        return indices;
      }
      indices[Dim] = 0;
      return do_next(indices, int_c<Dim + 1>);
    }

    constexpr optional<std::array<index, rank>>
    do_next(std::array<index, rank> indices, int_constant<rank - 1>) const
        noexcept {
      ++indices[rank - 1];
      if (indices[rank - 1] < extents().get(rank - 1)) {
        return indices;
      }
      return {};
    }

    template <int... Is>
    constexpr index do_mapping(const std::array<index, rank>& is,
                               std::integer_sequence<int, Is...>) const
        noexcept {
      std::array<index, rank> strides{{stride<Is>()...}};
      index sum = fub::transform_reduce(strides, is, index(0));
      return sum;
    }

    template <int... Rs>
    constexpr index do_stride(std::integer_sequence<int, Rs...>) const
        noexcept {
      std::array<index, sizeof...(Rs)> exts{{extents().get(Rs)...}};
      index prod = fub::accumulate(exts, index(1), std::multiplies<>());
      return prod;
    }
  };
};

} // namespace fub

#endif
