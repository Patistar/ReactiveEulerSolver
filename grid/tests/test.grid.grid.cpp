// Copyright (c) 2018 Maikel Nadolski
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "fub/serial/grid.hpp"
#include <type_traits>

#include "fub/optional.hpp"

struct Density {};

struct Equation {
  using complete_state = std::tuple<Density>;
};

int main()
{
    using fub::ready_future;
    using fub::serial::grid_node;
    using fub::extents;
    using fub::dyn;
    static_assert(std::is_assignable<grid_node<Equation, extents<dyn>>&, const grid_node<Equation, extents<dyn>>&>::value, "grid_node is not copyable.");

    ready_future<grid_node<Equation, extents<dyn>>> future_node(grid_node<Equation, extents<dyn>>(fub::serial::dummy_location{}, extents<dyn>(32)));
    grid_node<Equation, extents<dyn>> node = std::move(future_node);
    assert(node.get_patch_view().get().extents().size() == 32);
}
