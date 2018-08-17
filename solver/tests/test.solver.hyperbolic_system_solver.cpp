#include "fub/hyperbolic_system_solver.hpp"

#include "fub/euler/boundary_condition/reflective.hpp"
#include "fub/euler/hlle_riemann_solver.hpp"
#include "fub/euler/perfect_gas.hpp"
#include "fub/godunov_method.hpp"
#include "fub/output/cgns.hpp"
#include "fub/p4est/grid.hpp"
#include "fub/time_integrator/forward_euler.hpp"

#include <mpi.h>

#include <chrono>

static constexpr int Rank = 2;

using Equation = fub::euler::perfect_gas<Rank>;
using VariableList = fub::euler::perfect_gas_complete<Rank>;
using Grid = fub::p4est::basic_grid<VariableList, Rank, double>;
using Extents = fub::dynamic_extents_t<Rank>;

Grid make_simple_grid(Extents extents,
                      const fub::p4est::connectivity<Rank>& conn) {
  fub::uniform_cartesian_coordinates<Rank> coordinates({0.0, 0.0}, {1.0, 1.0},
                                                       fub::as_array(extents));
  fub::p4est::forest<Rank> forest(MPI_COMM_WORLD, conn, 0, 0, 1);
  fub::p4est::ghost_layer<Rank> ghost_layer(forest);
  return Grid(std::move(forest), std::move(ghost_layer),
              {VariableList(), extents, coordinates});
}

void fill_initial_values(Equation eq, fub::patch_t<Equation> data,
                         fub::coordinates_t<Grid> coords) {
  using namespace fub;
  using namespace fub::euler;
  for_each_index(data.get_mapping(), [=](std::ptrdiff_t i, std::ptrdiff_t j) {
    std::array<double, 2> x = coords.cell(i, j);
    x[0] -= 0.5;
    x[1] -= 0.5;
    if (fub::transform_reduce(x, x, 0.0) < 0.05) {
      data[density](i, j) = 1.;
      data[momentum<Rank, 0>](i, j) = 0.;
      data[momentum<Rank, 1>](i, j) = 0.;
      data[pressure](i, j) = 8.;
    } else {
      data[density](i, j) = 0.125;
      data[momentum<Rank, 0>](i, j) = 0.;
      data[momentum<Rank, 1>](i, j) = 0.;
      data[pressure](i, j) = 1.;
    }
    eq.reconstruct_prim(data(i, j), data(i, j));
  });
}

double L2_derivatives(int axis, Grid::const_patch_view left,
                      Grid::const_patch_view middle,
                      Grid::const_patch_view right) {
  using fub::euler::density;
  double sum = 0;
  for_each_index(middle.get_mapping(), [&](std::ptrdiff_t i, std::ptrdiff_t j) {
    const std::array<std::ptrdiff_t, 2> n{i, j};
    if (n[axis] == 0) {
      if (left.span()) {
        std::ptrdiff_t last_i = left.get_extents().extent(axis) - 1;
        std::array<std::ptrdiff_t, 2> last{n};
        last[axis] = last_i;
        sum += std::abs(fub::apply(middle[density], n) -
                        fub::apply(left[density], last));
        std::array<std::ptrdiff_t, 2> plast{last};
        plast[axis] -= 1;
        sum += std::abs(fub::apply(left[density], last) -
                        fub::apply(left[density], plast));
      }
      if (right.span()) {
        std::ptrdiff_t last_i = middle.get_extents().extent(axis) - 1;
        std::array<std::ptrdiff_t, 2> last{n};
        last[axis] = last_i;
        std::array<std::ptrdiff_t, 2> snd{n};
        snd[axis] = 1;
        sum += std::abs(fub::apply(middle[density], last) -
                        fub::apply(right[density], n));
        sum += std::abs(fub::apply(right[density], snd) -
                        fub::apply(right[density], n));
      }
    } else {
      sum += std::abs(middle[density](i, j) - middle[density](i - 1, j));
    }
  });
  std::ptrdiff_t volume = size(middle.get_extents());
  return sum / volume;
}

double L2_derivative(const Grid& grid, fub::quadrant_t<Grid> quad,
                     fub::axis axis) {
  fub::patch_buffer_t<Grid> left_buffer(grid.patch_extents());
  fub::patch_buffer_t<Grid> right_buffer(grid.patch_extents());
  const fub::face fl{axis, fub::direction::left};
  const fub::face fr{axis, fub::direction::right};
  Grid::const_patch_view middle = grid.patch_data(quad);
  Grid::const_patch_view left =
      interpolate_neighbor_data(grid, quad, fl, left_buffer.span());
  Grid::const_patch_view right =
      interpolate_neighbor_data(grid, quad, fr, right_buffer.span());
  return L2_derivatives(as_int(axis), left, middle, right);
}

bool is_large_L2_derivative(const Grid& grid, fub::quadrant_t<Grid> quad,
                            fub::axis axis) {
  const std::array<double, 2> dx = grid.coordinates(quad).dx();
  const double volume = fub::accumulate(dx, 1.0, std::multiplies<>{});
  return L2_derivative(grid, quad, axis) > 0.1 * volume;
}

bool is_large_L2_derivative(const Grid& grid, fub::quadrant_t<Grid> quad) {
  return is_large_L2_derivative(grid, quad, fub::axis::x) ||
         is_large_L2_derivative(grid, quad, fub::axis::y);
}

bool is_small_L2_derivative(const Grid& grid, fub::quadrant_t<Grid> quad,
                            fub::axis axis) {
  const std::array<double, 2> dx = grid.coordinates(quad).dx();
  const double volume = fub::accumulate(dx, 1.0, std::multiplies<>{});
  return L2_derivative(grid, quad, axis) < 0.01 * volume;
}

bool is_small_L2_derivative(const Grid& grid, fub::quadrant_t<Grid> quad) {
  return is_small_L2_derivative(grid, quad, fub::axis::x) &&
         is_small_L2_derivative(grid, quad, fub::axis::y);
}

void fill_initial_values(Equation eq, Grid& grid, int max_level) {
  using namespace fub;

  for_each_patch_box(grid, [&](auto&& box) {
    fill_initial_values(eq, grid.patch_data(box), grid.coordinates(box));
  });

  for (int i = 0; i < max_level; ++i) {
    grid.exchange_ghost_data();

    fub::p4est::forest<2> forest = grid.get_forest();

    refine_if(forest, fub::p4est::max_level(max_level),
              [&](int which_tree, fub::p4est::quadrant<2> quad) {
                return is_large_L2_derivative(grid, quad, axis::x) ||
                       is_large_L2_derivative(grid, quad, axis::y);
              });

    balance(forest);

    grid.interpolate_data_to(forest, [&](auto /* out */, auto in) {
      int size = in.quadrants.size();
      for (int i = 0; i < size; ++i) {
        fill_initial_values(eq, in.datas[i], grid.coordinates(in.quadrants[i]));
      }
    });

    partition(forest);
    grid.transfer_data_to(forest);
  }
}

struct State {
  Grid grid;
  std::chrono::duration<double> time;
  std::ptrdiff_t cycle;
};

void write_output(const State& state) {
  try {
    fub::output::cgns writer{};
    std::string filename =
        fmt::format("test.solver.hyperbolic_system_solver.{}.{}.cgns",
                    state.grid.get_forest().mpi_rank(), state.cycle);
    auto file = writer.open(filename.c_str(), Rank);
    writer.write(file, state);
  } catch (fub::output::cgns_error& error) {
    fmt::format("[CGNS] Error: {}\n", error.what());
  }
}

void refine(Grid& grid, int max_level) {
  grid.exchange_ghost_data();
  fub::p4est::forest<2> forest = grid.get_forest();

  refine_if(forest, fub::p4est::max_level(max_level),
            [&](int which_tree, fub::p4est::quadrant<2> quad) {
              return is_large_L2_derivative(grid, quad, fub::axis::x) ||
                     is_large_L2_derivative(grid, quad, fub::axis::y);
            });
  balance(forest);
  grid.interpolate_data_to(forest, [&](auto out, auto in) {
    assert(out.quadrants.size() <= in.quadrants.size());
    if (out.quadrants.size() < in.quadrants.size()) {
      assert(out.quadrants.size() == 1);
      for (int i = 0; i < in.quadrants.size(); ++i) {
        fub::linearily_refine(out.datas[0], in.datas[i],
                              child_id(in.quadrants[i]));
      }
    } else {
      assert(out.quadrants.size() == 1 && in.quadrants.size() == 1);
      std::memcpy(in.datas[0].span().data(), out.datas[0].span().data(),
                  out.datas[0].span().byte_size());
    }
  });

  grid.exchange_ghost_data();
  coarsen_if(forest,
             [&](int which_tree, fub::span<const fub::quadrant_t<Grid>> fine) {
               for (fub::quadrant_t<Grid> quad : fine) {
                 if (!(is_small_L2_derivative(grid, quad, fub::axis::x) &&
                       is_small_L2_derivative(grid, quad, fub::axis::y))) {
                   return false;
                 }
               }
               return true;
             });
  balance(forest);

  grid.interpolate_data_to(forest, [&](auto out, auto in) {
    assert(out.quadrants.size() >= in.quadrants.size());
    if (out.quadrants.size() > in.quadrants.size()) {
      assert(in.quadrants.size() == 1);
      for (int i = 0; i < out.quadrants.size(); ++i) {
        fub::linearily_coarsen(out.datas[i], in.datas[0],
                               child_id(out.quadrants[i]));
      }
    } else {
      assert(out.quadrants.size() == 1 && in.quadrants.size() == 1);
      std::memcpy(in.datas[0].span().data(), out.datas[0].span().data(),
                  out.datas[0].span().byte_size());
    }
  });

  partition(forest);
  grid.transfer_data_to(forest);
}

void test_solver() {
  using namespace fub;

  // Initialize Application

  fub::p4est::connectivity<Rank> conn = fub::p4est::unit_square;
  const int max_level = 7;
  State state{make_simple_grid(Extents(8, 8), conn)};

  auto solver = make_hyperbolic_system_solver(
      godunov_method<euler::hlle_riemann_solver>(),
      time_integrator::forward_euler(), euler::perfect_gas<Rank>(1.4, 28.),
      euler::boundary_condition::reflective());

  // Reapply Initial Conditions and Refine Grid

  fill_initial_values(solver.equation, state.grid, max_level);

  write_output(state);

  // Do Time Stepping

  int n_cycles = 1000;
  std::chrono::duration<double> final_time(1.0);
  std::chrono::duration<double> eps(std::numeric_limits<double>::epsilon());

  while (state.time < final_time) {
    if (state.cycle > 0 && state.cycle % 2 == 0) {
      refine(state.grid, max_level);
    }

    const std::chrono::duration<double> step_size_limit =
        final_time - state.time + 4 * eps;

    const std::chrono::duration<double> stable_time_step_size = [&] {
      const double dt =
          solver.estimate_stable_time_step_size(state.grid).count();
      double count = std::numeric_limits<double>::signaling_NaN();
      MPI_Allreduce(&dt, &count, 1, MPI_DOUBLE, MPI_MIN,
                    state.grid.get_forest().mpi_communicator());
      return std::chrono::duration<double>(count);
    }();

    const std::chrono::duration<double> time_step_size =
        std::min(stable_time_step_size, step_size_limit);

    state.grid = solver.time_step(state.grid, time_step_size);
    state.time += time_step_size;
    state.cycle += 1;

    if (state.grid.get_forest().mpi_rank() == 0) {
      int percentage = std::ceil((100 * state.time) / final_time);
      fmt::print("[{:3}\%] time: {}s time_step_size: {}s\n", percentage,
                 state.time.count(), time_step_size.count());
    }

    write_output(state);
  }
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  try {
    test_solver();
  } catch (...) {
  }
  MPI_Finalize();
}