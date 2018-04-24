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

#include "fub/euler/ideal_gas.hpp"
#include "fub/euler/mechanism/single_stage.hpp"
#include "fub/euler/variables.hpp"
#include "fub/output/cgns.hpp"
#include "fub/patch.hpp"
#include "fub/span.hpp"
#include <iostream>

using fub::euler::Density;
using fub::euler::ideal_gas;
using fub::euler::mechanism::single_stage::single_stage;
using fub::output::cgns;

int main() {
  fub::patch<ideal_gas<single_stage, 2>::complete_state, fub::extents<32, 32>>
      patch{};
  std::fill(patch.get<Density>().begin(),
            patch.get<Density>().begin() + 32 * 16, 0.125);
  std::fill(patch.get<Density>().begin() + 32 * 16, patch.get<Density>().end(),
            1.0);
  fub::uniform_cartesian_coordinates<2> coordinates{
      {0.0, 0.0}, {1, 1}, {32, 32}};
  try {
    auto ctx = cgns::open("test.cgns", 2);
    cgns::write(ctx, fub::octant<2>{1, {0, 0}}, fub::make_view(patch),
                coordinates, ideal_gas<single_stage, 2>());
    cgns::write(ctx, fub::octant<2>{1, {1, 0}}, fub::make_view(patch),
                coordinates, ideal_gas<single_stage, 2>());
    cgns::write(ctx, fub::octant<2>{1, {0, 1}}, fub::make_view(patch),
                coordinates, ideal_gas<single_stage, 2>());
    cgns::write(ctx, fub::octant<2>{1, {1, 1}}, fub::make_view(patch),
                coordinates, ideal_gas<single_stage, 2>());
  } catch (std::exception& e) {
    std::cerr << e.what() << '\n';
  }
}
