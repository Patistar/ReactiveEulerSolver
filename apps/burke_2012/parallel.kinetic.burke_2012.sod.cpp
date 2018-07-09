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

#include "fub/parallel/kinetic.burke_2012.1d.hpp"

#include "fub/euler/boundary_condition/reflective.hpp"
#include "fub/output/cgns.hpp"
#include "fub/patch_view.hpp"
#include "fub/run_simulation.hpp"
#include "fub/uniform_cartesian_coordinates.hpp"

#include <boost/program_options.hpp>

#include <array>

using Equation = fub::parallel::kinetic::burke_2012_1d::equation_type;
using Grid = fub::parallel::kinetic::burke_2012_1d::grid_type;
using Partition = Grid::partition_type;

std::array<Equation::complete_state, 2> get_initial_states() noexcept {
  std::array<double, Equation::species_size> moles{};
  using namespace fub::euler::mechanism::burke2012::variables;
  moles[as_index(n2)] = 1.0;
  auto left = Equation().set_TPX(300, 10E5, moles);
  std::swap(moles[as_index(n2)], moles[as_index(ar)]);
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

using state_type = fub::parallel::kinetic::burke_2012_1d::state_type;

struct write_cgns_file {
  mutable hpx::future<void> queue = hpx::make_ready_future();

  void operator()(state_type state) const {
    fmt::print("[CGNS] Add cgns output to the queue...\n", state.cycle);
    queue = queue.then([state = std::move(state)](hpx::future<void>) {
      std::string file_name = fmt::format("out_{}.cgns", state.cycle);
      auto file = fub::output::cgns::open(file_name.c_str(), 2);
      fub::output::cgns::iteration_data_write(file, state.time, state.cycle);
      for (const Partition& partition : state.grid) {
        const auto& octant = fub::grid_traits<Grid>::octant(partition);
        auto data = partition.second.get_patch_view().get();
        fub::output::cgns::write(file, octant, data, state.coordinates,
                                 Equation());
      }
    });
  }

  ~write_cgns_file() noexcept {
    fmt::print("[CGNS] Wait for all output tasks to finish their job.\n");
    queue.wait();
  }
};

int main(int argc, char** argv) {
  namespace po = boost::program_options;
  po::options_description desc("Allowed Options");
  desc.add_options()("depth", po::value<int>()->default_value(3),
                     "Depth of tree.");
  desc.add_options()("time", po::value<double>()->default_value(1.0),
                     "The final time level which we are interested in.");
  desc.add_options()("feedback_interval",
                     po::value<double>()->default_value(1e-4),
                     "The time interval in which we write output files.");
  desc.add_options()("extents", po::value<fub::index>()->default_value(12),
                     "Amount of cells per patch.");
  return hpx::init(desc, argc, argv);
}

int hpx_main(boost::program_options::variables_map& vm) {
  const int depth = vm["depth"].as<int>();
  
  const std::array<fub::index, 1> extents{{vm["extents"].as<fub::index>()}};
  fub::uniform_cartesian_coordinates<1> coordinates({0}, {1.0}, extents);
  
  auto state = fub::parallel::kinetic::burke_2012_1d::initialise(
      &initial_value_function, coordinates, depth);
  write_cgns_file write_cgns{};
  write_cgns(state);
  fub::euler::boundary_condition::reflective boundary_condition{};
  fub::simulation_options options{};
  options.feedback_interval =
      std::chrono::duration<double>(vm["feedback_interval"].as<double>());
  options.final_time = std::chrono::duration<double>(vm["time"].as<double>());
  fub::run_simulation(fub::parallel::kinetic::burke_2012_1d(), state,
                      boundary_condition, options, std::ref(write_cgns),
                      fub::print_cycle_timings{options.final_time});
  return hpx::finalize();
}
