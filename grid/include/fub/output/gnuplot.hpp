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

#ifndef FUB_OUTPUT_GNUPLOT_HPP
#define FUB_OUTPUT_GNUPLOT_HPP

#include "fub/tuple.hpp"
#include <fmt/format.h>
#include <array>

namespace fub {
inline namespace v1 {
namespace output {

struct gnuplot {
  template <typename Patch, typename CoordinateMapping>
  static void write(std::FILE* file, const Patch& patch, CoordinateMapping x) {
    fmt::print(file, "#");
    for (int dim = 0; dim < x.rank(); ++dim) {
      fmt::basic_memory_buffer<char, 4> buffer;
      fmt::format_to(buffer, "X({})", dim);
      std::array<char, 5> array{};
      std::copy_n(buffer.data(), 4, array.data());
      fmt::print(file, "{:>15}", array.data());
    }
    for_each(patch.get_variable_list(), [&](auto variable) {
      fmt::print(file, "{:>15}", variable.name());
    });
    fmt::print(file, "\n ");
    for_each_index(patch.get_mapping(), [&](auto... is) {
      auto xs = x.cell(is...);
      for (auto x_i : xs) {
        fmt::print(file, "{:>15.6f}", x_i);
      }
      for_each(patch.get_variable_list(), [&](auto variable) {
        fmt::print(file, "{:>#15.10f}", patch[variable](is...));
      });
      fmt::print(file, "\n ");
    });
    fmt::print(file, "\n");
  }

  template <typename State>
  static void write(std::FILE* file, const State& state) {
    fmt::print(file, "# time = {:<.6e}, cycle = {}\n", state.time.count(),
               state.cycle);
    auto local_nodes = state.grid.get_local_nodes();
    for (auto&& node : local_nodes) {
      gnuplot::write(file, state.grid.get_patch_data(node),
                     state.grid.get_coordinates(node));
    }
  }
};

} // namespace output
} // namespace v1
} // namespace fub

#endif // !FUB_OUTPUT_GNUPLOT_HPP
