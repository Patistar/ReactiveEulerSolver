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

#ifndef FUB_PARALLEL_INITIALISE_HPP
#define FUB_PARALLEL_INITIALISE_HPP

#include "fub/parallel/grid.hpp"
#include "fub/uniform_cartesian_coordinates.hpp"

namespace fub {
namespace parallel {
template <typename Equation, typename Extents>
grid<Equation, Extents> make_grid(int depth, const Equation& equation,
                                  const Extents& extents) {
  grid<Equation, Extents> grid(equation, extents);
  octant<Extents::rank> oct{depth, {}};
  octant<Extents::rank> last = octant<Extents::rank>{}.upper_descendant(depth);
  while (oct != last) {
    grid.insert(grid.end(), oct,
                grid_node<Equation, Extents>(dummy_location{}, extents));
    oct = oct.next();
  }
  grid.insert(grid.end(), oct,
              grid_node<Equation, Extents>(dummy_location{}, extents));
  return grid;
}

template <typename State, typename InitialCondition, int Rank>
State initialise(InitialCondition initial_condition,
                 const uniform_cartesian_coordinates<Rank>& coordinates,
                 int depth, double cfl = 0.5) {
  using grid_type = std::decay_t<decltype(std::declval<State>().grid)>;
  using traits = grid_traits<grid_type>;
  typename traits::extents_type extents(coordinates.extents());
  typename traits::equation_type equation{};
  grid_type grid = make_grid(depth, equation, extents);
  for (typename traits::partition_type& partition : grid) {
    auto&& octant = traits::octant(partition);
    auto&& coords = adapt(coordinates, octant);
    auto&& node = traits::node(partition);
    node = traits::dataflow_action(
        [](typename traits::node_type node,
           const uniform_cartesian_coordinates<Rank>& coords,
           InitialCondition initial_condition) {
          return node.get_patch_view().then([node, coords,
                                             initial_condition](auto view) {
            auto patch_view = view.get();
            typename traits::patch_type patch(patch_view.extents());
            const auto pv = make_view(patch);
            for_each_index(pv.extents(), [&](const std::array<index, Rank>& i) {
              pv(i) = fub::invoke(initial_condition, fub::apply(coords, i));
            });
            return
                typename traits::node_type(dummy_location{}, std::move(patch));
          });
        },
        dummy_location(), node, coords, initial_condition);
  }
  return State{std::move(grid),
               coordinates,
               std::chrono::duration<double>(0),
               std::chrono::duration<double>(0),
               0,
               cfl};
}

} // namespace parallel
} // namespace fub

#endif // !FUB_PARALLEL_INITIALISE_HPP
