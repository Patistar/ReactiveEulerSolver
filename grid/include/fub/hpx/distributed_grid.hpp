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

#ifndef DISTRIBUTED_GRID_HPP
#define DISTRIBUTED_GRID_HPP

#if defined(__CLANG__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#endif
#include <hpx/include/components.hpp>
#if defined(__CLANG__)
#pragma clang diagnostic pop
#endif

#include "fub/interval_map.hpp"
#include "fub/optional.hpp"

namespace fub {
namespace hpx {

template <typename Eq, typename E>
class remote_grid_client
    : ::hpx::components::client_base<remote_grid_client, remote_grid> {};

template <typename Eq, typename E> class distributed_grid {
public:
  using client_type = remote_grid_client<Eq, E>;

  explicit distributed_grid(const interval_map<octant<rank>, hpx::id_type>&);

  client_type* find_remote_client(octant<rank>);

private:
  interval_map<octant<rank>, client_type> m_remote_grids;
};

} // namespace hpx
} // namespace fub

#endif // !DISTRIBUTED_GRID_HPP
