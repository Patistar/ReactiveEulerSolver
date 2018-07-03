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

#ifndef FUB_DRIVER_DISTRIBUTED_INITIALISE_HPP
#define FUB_DRIVER_DISTRIBUTED_INITIALISE_HPP

#include "fub/distributed/grid.hpp"
#include "fub/uniform_cartesian_coordinates.hpp"

namespace fub {
namespace distributed {
template <typename Grid, typename Function> struct initialise_patch {
  using traits = grid_traits<Grid>;
  using node_type = typename traits::node_type;
  using patch_type = typename traits::patch_type;
  using extents_type = typename traits::extents_type;
  static constexpr int rank = traits::rank;

  static node_type invoke(extents_type extents,
                          uniform_cartesian_coordinates<rank> coords,
                          Function function) {
    patch_type patch(extents);
    auto view = make_view(patch);
    for_each_index(view.extents(), [&](const std::array<index, rank>& i) {
      view(i) = fub::invoke(function, fub::apply(coords, i));
    });
    return node_type(hpx::find_here(), std::move(patch));
  }
};

template <typename Grid, typename Feedback>
struct initialise_patch_action
    : hpx::actions::make_action<
          decltype(&initialise_patch<Grid, Feedback>::invoke),
          &initialise_patch<Grid, Feedback>::invoke,
          initialise_patch_action<Grid, Feedback>> {};

template <typename State, typename InitialCondition, int Rank>
State initialise(InitialCondition f,
                 const uniform_cartesian_coordinates<Rank>& coordinates,
                 int depth, double cfl = 0.5) {
  using grid_type = decltype(std::declval<State>().grid);
  using traits = grid_traits<grid_type>;
  typename traits::extents_type extents(coordinates.extents());
  typename traits::equation_type equation{};
  grid_type grid = make_grid(make_uniform_distribution(int_c<Rank>, depth, hpx::find_all_localities()), depth, extents, equation);
  for (auto& partition : grid) {
    hpx::id_type locality = traits::locality(partition);
    auto&& oct = traits::octant(partition);
    auto&& coords = adapt(coordinates, oct);
    partition.second =
        traits::dataflow_action(initialise_patch_action<grid_type, InitialCondition>(),
                   locality, extents, coords, f);
  }
  return State{std::move(grid),
               coordinates,
               std::chrono::duration<double>(0),
               std::chrono::duration<double>(0),
               0,
               cfl};
}

} // namespace distributed
} // namespace fub

#endif // !FUB_PARALLEL_INITIALISE_HPP
