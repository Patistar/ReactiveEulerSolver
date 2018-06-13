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

#ifndef INITIALISE_HPP
#define INITIALISE_HPP

#include "fub/functional.hpp"
#include "fub/grid.hpp"
#include "fub/uniform_cartesian_coordinates.hpp"

namespace fub {
namespace serial {

template <typename State, typename InitialCondition, int Rank>
State initialise(InitialCondition f,
                 const uniform_cartesian_coordinates<Rank>& coordinates,
                 int depth, double cfl = 0.5) {
  using grid_type = decltype(std::declval<State>().grid);
  using traits = grid_traits<grid_type>;
  typename grid_type::extents_type extents(coordinates.extents());
  grid_type grid(depth, extents);
  for (auto& partition : grid) {
    auto octant = traits::octant(partition);
    partition.second = traits::dataflow(
        [=](typename traits::patch_type patch) {
          auto adapted = adapt(coordinates, octant);
          auto view = make_view(patch);
          for_each_index(patch.extents(),
                         [&](const std::array<index, Rank>& i) {
                           view(i) = fub::invoke(f, fub::apply(adapted, i));
                         });
          return std::make_shared<typename traits::node_type>(std::move(patch));
        },
        std::move(partition));
  }
  return State{std::move(grid),
               coordinates,
               std::chrono::duration<double>(0),
               std::chrono::duration<double>(0),
               0,
               cfl};
}

} // namespace serial
} // namespace fub

#endif // !INITIALISE_HPP
