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

#include "fub/serial/kinetic.gri_30.1d.hpp"

#include "fub/euler/boundary_condition/reflective.hpp"
#include "fub/grid.hpp"
#include "fub/output/cgns.hpp"
#include "fub/patch_view.hpp"
#include "fub/run_simulation.hpp"
#include "fub/uniform_cartesian_coordinates.hpp"

#include <boost/program_options.hpp>

#include <array>

using Equation = fub::serial::kinetic::gri_30_1d::equation_type;
using Grid = fub::serial::kinetic::gri_30_1d::grid_type;
using Partition = Grid::partition_type;

std::array<Equation::complete_state, 2> get_initial_states() noexcept {
  std::array<double, Equation::species_size> moles{};
  using namespace fub::euler::mechanism::gri_30;
  moles[as_index(CH4())] = 1.0;
  moles[as_index(O2())] = 1.0;
  moles[as_index(C2H6())] = 1.0;
  auto left = Equation().set_TPX(600, 10E5, moles);
  std::swap(moles[as_index(N2())], moles[as_index(AR())]);
  auto right = Equation().set_TPX(345, 1E5, moles);
  return {{left, right}};
}

Equation::complete_state
initial_value_function(const std::array<double, 1>& xs) {
  static auto states = get_initial_states();
  if (xs[0] < 0.1) {
    return states[0];
  } else {
    return states[1];
  }
}

using state_type = fub::serial::kinetic::gri_30_1d::state_type;

struct write_cgns_file {
  void operator()(state_type state) const {
    fmt::print("[CGNS] Write output file...\n", state.cycle);
    std::string file_name = fmt::format("out_{}.cgns", state.cycle);
    auto file = fub::output::cgns::open(file_name.c_str(), 2);
    fub::output::cgns::iteration_data_write(file, state.time, state.cycle);
    for (const Partition& partition : state.grid) {
      const auto& octant = fub::grid_traits<Grid>::octant(partition);
      auto node = partition.second.get();
      fub::output::cgns::write(file, octant, fub::make_view(node->patch),
                               state.coordinates, Equation());
    }
  }
};

int main(int argc, char** argv) {
  namespace po = boost::program_options;
  po::options_description desc("Allowed Options");
  desc.add_options()("depth", po::value<int>()->default_value(0),
                     "Depth of tree.");
  desc.add_options()("time", po::value<double>()->default_value(1e-4),
                     "The final time level which we are interested in.");
  desc.add_options()("feedback_interval",
                     po::value<double>()->default_value(1e-6),
                     "The time interval in which we write output files.");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  const int depth = vm["depth"].as<int>();
  auto extents = static_cast<fub::array<fub::index, 1>>(Grid::extents_type());

  fub::uniform_cartesian_coordinates<1> coordinates({0}, {1.0}, extents);

  auto state = fub::serial::kinetic::gri_30_1d::initialise(
      &initial_value_function, coordinates, depth);
  write_cgns_file write_cgns;
  write_cgns(state);
  fub::euler::boundary_condition::reflective boundary_condition{};
  fub::simulation_options options{};
  options.feedback_interval =
      std::chrono::duration<double>(vm["feedback_interval"].as<double>());
  options.final_time = std::chrono::duration<double>(vm["time"].as<double>());
  fub::run_simulation(fub::serial::kinetic::gri_30_1d(), state,
                      boundary_condition, options, write_cgns,
                      fub::print_cycle_timings{options.final_time});
}
