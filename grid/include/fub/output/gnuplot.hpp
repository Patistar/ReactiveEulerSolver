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

#ifndef FUB_GNUPLOT_HPP
#define FUB_GNUPLOT_HPP

#include "fub/algorithm.hpp"
#include "fub/octree.hpp"
#include "fub/patch_view.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <functional>
#include <iomanip>
#include <ostream>

namespace fub {
namespace output {

struct Gnuplot {
  template <typename I, std::size_t Rank>
  static void write(std::ostream& out, const std::array<I, Rank>& array) {
    out << '{';
    out << array[0];
    for (std::size_t i = 1; i < Rank; ++i) {
      out << ", " << array[i];
    }
    out << '}';
  }

  template <int Rank, typename CoordinateMapping, typename Extent,
            typename... Quantities>
  static void write(std::ostream& out, const octant<Rank>& octant,
                    const View<Extent, Quantities...>& partition,
                    CoordinateMapping coords, std::size_t cycle,
                    std::chrono::duration<double> t) {
    out << "# Time = " << t.count() << ", Cycle = " << cycle << ", octant: ";
    write(out, coordinates(octant));
    std::bitset<64> morton(octant.morton_index());
    out << " (" << morton << ")\n# ";
    for (int dim = 0; dim < Rank; ++dim) {
      out << std::setw(11) << "X(" << dim << ") ";
    }
    fub::for_each_tuple_element(
        [&out](auto quantity) {
          out << std::setw(15) << quantity.name() << ' ';
        },
        std::tuple<Quantities...>());
    out << '\n';
    fub::for_each_index(
        partition.extents(), [&](const std::array<index, Rank>& is) {
          auto xs = fub::apply(coords, is);
          for (auto x : xs) {
            out << std::setw(15) << x << ' ';
          }
          fub::for_each_tuple_element(
              [&](auto quantity) {
                out << std::setw(15) << fub::apply(partition[quantity], is)
                    << ' ';
              },
              std::tuple<Quantities...>());
          out << '\n';
        });
    out << '\n';
  }
};

} // namespace output
} // namespace fub

#endif // !GNUPLOT_HPP
