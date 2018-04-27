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

#ifndef FUB_DRIVER_SERIAL_GRI_30_1D_HPP
#define FUB_DRIVER_SERIAL_GRI_30_1D_HPP

#include "fub/grid.hpp"

#include "fub/euler/mechanism/gri_30.hpp"
#include "fub/polymorphic_boundary_condition.hpp"
#include "fub/uniform_cartesian_coordinates.hpp"

#include <chrono>
#include <cstdint>
#include <functional>

namespace fub {
namespace serial {
namespace kinetic {

struct gri_30_1d {
  static constexpr int rank = 1;
  static constexpr int ghost_width = 2;

  using patch_extents_type = extents<200>;
  using equation_type =
      euler::ideal_gas<euler::mechanism::gri_30::gri_30, rank>;
  using grid_type = standard_grid<equation_type, patch_extents_type>;
  using partition_type = grid_type::partition_type;
  using patch_type = grid_type::patch_type;
  using equation_state_t = complete_state_t<equation_type>;

  struct state_type {
    grid_type grid;
    uniform_cartesian_coordinates<rank> coordinates;
    std::chrono::duration<double> time;
    std::chrono::duration<double> dt;
    std::ptrdiff_t cycle;
    double cfl;
  };

  using initial_condition_function =
      std::function<equation_state_t(const std::array<double, rank>&)>;

  using feedback_function = std::function<bool(const state_type&)>;

  using boundary_condition =
      polymorphic_boundary_condition<grid_type, ghost_width, axis::x>;

  static state_type initialise(initial_condition_function,
                               const uniform_cartesian_coordinates<rank>&,
                               int depth);

  static state_type advance(const state_type& state,
                            std::chrono::duration<double> dt,
                            const boundary_condition& boundary,
                            feedback_function feedback);
};

} // namespace kinetic
} // namespace serial
} // namespace fub

#endif
