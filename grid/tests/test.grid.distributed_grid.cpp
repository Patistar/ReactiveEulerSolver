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

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#endif
#include <parallel/hpx.hpp>
#ifdef __clang__
#pragma clang diagnostic pop
#endif

#include "fub/distributed/grid.hpp"

#include "fub/extents.hpp"
#include "fub/interval_map.hpp"
#include "fub/variables.hpp"

#include <cassert>
#include <fmt/format.h>

struct Density {};

struct Equation {
  using complete_state = std::tuple<Density>;
};

using Extents = ::fub::extents<fub::dyn>;

FUB_REGISTER_GRID_NODE_COMPONENT(Equation, Extents);

void print_distribution(
    const fub::interval_map<fub::octant<Extents::rank>, hpx::id_type>& dist) {
  const int depth = 3;
  fub::octant<Extents::rank> o{depth, {0}};
  const fub::octant<Extents::rank> last =
      fub::octant<Extents::rank>().upper_descendant(depth);
  while (o != last) {
    std::bitset<64> morton(o.morton_index());
    fmt::print("Ocant: {} -> hpx::id_type: {}\n", morton, dist[o]);
    o = o.next();
  }
  std::bitset<64> morton(o.morton_index());
  fmt::print("Ocant: {} -> hpx::id_type: {}\n", morton, dist[o]);
}

int hpx_main() {
  // Make a distribution
  static constexpr int rank = Extents::rank;
  const int depth = 3;
  const auto localities = hpx::find_all_localities();
  const auto distribution = fub::distributed::make_uniform_distribution(
      fub::int_c<rank>, depth, localities);

  print_distribution(distribution);

  // Make a grid
  Extents extents{32};
  Equation equation{};
  const auto grid =
      fub::distributed::make_grid(distribution, depth, extents, equation);
  return hpx::finalize();
}

int main(int argc, char** argv) { hpx::init(argc, argv); }
