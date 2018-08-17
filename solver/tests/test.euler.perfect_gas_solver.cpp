// Copyright(c) 2018 Maikel Nadolski
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

#include "fub/euler/hlle_riemann_solver.hpp"
#include "fub/euler/perfect_gas.hpp"
#include "fub/godunov_method.hpp"
#include "fub/hyperbolic_system_solver.hpp"
#include "fub/p4est/grid.hpp"

#include <boost/program_options.hpp>

namespace po = boost::program_options;

constexpr int Rank = 2;

using equation_t = fub::euler::perfect_gas<Rank>;

using grid_t = fub::p4est::basic_grid<equation_t::complete_state, Rank>;

using simulation_state_t = fub::simulation_state<grid_t, equation_t>;

/// Returns a variables map which contains a run-time configuration.
po::variables_map Euler_parse_args(int argc, char** argv);

/// This is our real main function after some boot strapping like initialising
/// MPI and reading command line options.
void Euler_main(const po::variables_map& options);

/// Initialises a simulation state with specified program options
simulation_state_t Euler_initialise(const po::variables_map& options);

void Euler_write_output(const simulation_state_t& state);

void Euler_checkpoint(const simulation_state_t& state) noexcept;

/// Checks the error code and prints an MPI error message
void MPI_Abort_on_error(int error_code) noexcept;

int main(int argc, char** argv) {
  MPI_Abort_on_error(MPI_Init(&argc, &argv));
  try {
    Euler_main(Euler_parse_args(argc, argv));
  } catch (std::exception& error) {
    fmt::print(stderr, "An unexpected exception was thrown.\n");
    fmt::print(stderr, "error.what(): \"{}\"\n", error.what());
    fmt::print(stderr, "Terminating the program.\n")
  }
  MPI_Finalize();
}

void MPI_Abort_on_error(int error_code) noexcept {
  if (error_code != MPI_SUCCESS) {
    fmt::print(stderr, "An unexpected MPI error happened.\n");
    fmt::print(stderr, "what(): \"{}\"\n");
    fmt::print(stderr, "Terminating the program.\n");
    std::terminate();
  }
}

po::variables_map Euler_parse_args(int argc, char** argv) {
  po::options_description desc("Allowed Options");
  desc.add_options()("help", "produce help message")(
      "initial_level", po::value<int>(), "initial refinement level")(
      "max_level", po::value<int>(), "maximal refinement level")(
      "gamma", po::value<double>(), "heat capacity ratio of the equation")(
      "patch_extents", po::value<int>(),
      "amount of cells in each direction per patch")(
      "upper_x", po::value<double>(), "upper bound for x coordinates")(
      "upper_y", po::value<double>(), "upper bound for y coordinates");
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  return vm;
}

void Euler_main(po::variables_map options) {
  // Build our solver
  auto solver = fub::make_hyperbolic_system_solver(
      fub::godunov_method(fub::euler::hlle_riemann_solver()),
      fub::time_integrator::forward_euler());

  // Initialise the simulation state
  simulation_state_t state = Euler_initialise(options);

  // Get some simulation options
  double output_interval = options["output_interval"].as<double>();
  double max_time = options["max_time"].as<double>();
  double max_wall_time = options["max_wall_time"].as<double>();
  double max_cycle = options["max_cycles"].as<double>();
  fub::print_statistics_line cycle_feedback(stdout, max_time, max_cycle,
                                            max_wall_time);
  try {
    // Loop in time until we achieved our goal.
    while (state.time < max_time && state.cycle < max_cycle) {
      constexpr double eps = std::numeric_limits<double>::eps();
      const double dt = std::min(output_interval, max_time - state.time + eps);
      fub::regrid(state);
      fub::advance_time(state, solver, duration_t(dt), cycle_feedback);
      Euler_write_output(state);
    }
  } catch (...) {
    // An unsexpected exception happend while looping in time. Try to save the
    // last legal state and rethrow.
    Euler_checkpoint(state);
    throw;
  }
}

} // namespace example
} // namespace fub