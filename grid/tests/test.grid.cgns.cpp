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

#include "fub/output/cgns.hpp"

#include "fub/uniform_cartesian_coordinates.hpp"
#include "fub/variable_data.hpp"

#include <cstdio>

struct Density : fub::scalar_variable {
  static const char* name(int) noexcept { return "Density"; }
};
static constexpr auto density = fub::tag<Density>;

struct Momentum : fub::vector_variable<2> {
  static std::array<char, 10> construct_name(int i) noexcept {
    fmt::basic_memory_buffer<char, 10> buffer{};
    std::array<const char*, 3> names = {"X", "Y", "Z"};
    fmt::format_to(buffer, "Momentum{}", names[i % 3]);
    std::array<char, 10> array{};
    std::copy_n(buffer.data(), buffer.size(), array.begin());
    return array;
  }

  static auto construct_names(int n) noexcept {
    std::array<std::array<char, 10>, 2> arrays;
    for (int i = 0; i < n; ++i) {
      arrays[i] = construct_name(i);
    }
    return arrays;
  }

  static const char* name(int i) noexcept {
    static const std::array<std::array<char, 10>, 2> name_ = construct_names(2);
    return name_[i].data();
  }
};

using Variables = fub::variable_list<Density, Momentum>;

int main() {
  Variables vars{};
  fub::dynamic_extents_t<2> extents(1024, 1024);
  fub::variable_data<Variables, double, 2> patch(vars, extents);
  fub::uniform_cartesian_coordinates<2> coordinates({0., 0.}, {1., 1.},
                                                    fub::as_array(extents));
  for_each_index(patch.get_mapping(), [&](std::ptrdiff_t i, std::ptrdiff_t j) {
    auto x = coordinates.cell(i, j);
    if (x[0] + x[1] < 1.0) {
      patch[density](i, j) = 1.0;
    }
  });
  fub::output::cgns output_module;
  auto file = output_module.open("test.cgns", 2);
  output_module.write(file, patch, coordinates);
}
