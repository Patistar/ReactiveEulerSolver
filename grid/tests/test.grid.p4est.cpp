#include "fub/output/gnuplot.hpp"
#include "fub/p4est/grid.hpp"
#include <chrono>

using namespace fub;

struct Density : scalar_variable {
  static const char* name(int) noexcept { return "Density"; }
};
static constexpr auto density = tag<Density>;

struct Pressure : scalar_variable {
  static const char* name(int) noexcept { return "Pressure"; }
};
static constexpr auto pressure = tag<Pressure>;

struct Momentum : vector_variable<2> {
  static std::array<char, 10> construct_name(int i) noexcept {
    fmt::basic_memory_buffer<char, 10> buffer{};
    std::array<const char*, 3> names = {"X", "Y", "Z"};
    fmt::format_to(buffer, "Momentum{}", names[i % 3]);
    std::array<char, 10> array{};
    std::copy_n(buffer.data(), buffer.size(), array.begin());
    return array;
  }

  static auto construct_names(int n) noexcept {
    std::array<std::array<char, 10>, 2> arrays;
    for (int i = 0; i < n; ++i) {
      arrays[i] = construct_name(i);
    }
    return arrays;
  }

  static const char* name(int i) noexcept {
    static const std::array<std::array<char, 10>, 2> name_ = construct_names(2);
    return name_[i].data();
  }
};

using Variables = variable_list<Density, Momentum>;

struct State {
  v1::p4est::grid<Variables, 2> grid;
  std::chrono::duration<double> time;
  std::ptrdiff_t cycle;
};

void my_main(MPI_Comm communicator, int rank) {
  using namespace fub::v1::p4est;

  dynamic_extents_t<2> extents(4, 4);
  uniform_cartesian_coordinates<2> coordinates({0.0, 0.0}, {1.0, 1.0},
                                               as_array(extents));
  grid<Variables, 2> grid(communicator, unit_square, min_level(3), fill_uniform,
                          Variables(), extents, coordinates);

  auto trees = grid.get_local_trees();
  int counter = 0;
  for (auto& tree : trees) {
    for (auto& quad : tree.get_quadrants()) {
      auto coords = quad.get_coordinates();
      fmt::print("#{:<4}x: {}, y: {}, level: {}\n", rank, counter, coords[0],
                 coords[1], quad.get_level());
    }
    counter += 1;
  }
  fmt::print("#{:<4}Ghosts:\n", rank);
  for (auto& ghost : grid.get_ghost_quadrants()) {
    auto coords = ghost.get_coordinates();
    fmt::print("#{:<4}x: {}, y: {}, level: {}\n", rank, ghost.local_num(),
               ghost.which_tree(), coords[0], coords[1], ghost.get_level(),
               grid.get_owner(ghost));
  }
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  my_main(MPI_COMM_WORLD, rank);
  MPI_Finalize();
}