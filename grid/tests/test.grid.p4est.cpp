#include "fub/linearily_interpolate.hpp"
#include "fub/output/cgns.hpp"
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
template <int Dim> static constexpr auto momentum = tag<Momentum, Dim>;

using Variables = variable_list<Density, Momentum, Pressure>;

using Grid = fub::p4est::basic_grid<Variables, 2, double>;

void my_initialize(Grid::patch_view data, Grid::coordinates_type coords) {
  for_each_index(data.get_mapping(), [=](std::ptrdiff_t i, std::ptrdiff_t j) {
    std::array<double, 2> x = coords(i, j);
    if (x[0] + x[1] < 1.0) {
      data[density](i, j) = 1.;
      data[momentum<0>](i, j) = 0.;
      data[momentum<1>](i, j) = 0.;
      data[pressure](i, j) = 1.;
    } else {
      data[density](i, j) = 0.125;
      data[momentum<0>](i, j) = 0.;
      data[momentum<1>](i, j) = 0.;
      data[pressure](i, j) = 1.;
    }
  });
}

double L2_derivatives(int axis, Grid::const_patch_view left,
                      Grid::const_patch_view middle,
                      Grid::const_patch_view right) {
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

Grid::const_patch_view
interpolate_neighbor_data(const Grid& grid, fub::p4est::quadrant<2> quad,
                          fub::face face, fub::span<double> buffer) noexcept {

  span<const fub::p4est::mesh_quadrant<2>> neighbors =
      grid.face_neighbors(quad, face);

  // We hit a domain boundary
  if (neighbors.size() == 0) {
    return {};
  }
  int offset = static_cast<int>(face.side == direction::left);
  Grid::patch_view interpolated(Variables(), buffer, grid.patch_extents());

  if (neighbors[0].is_ghost) {
    return {};
  }

  // neighbors have a finer level than quad.
  if (neighbors.size() > 1) {
    assert(neighbors.size() == 2);
    assert(neighbors[0].quad.level() == quad.level() + 1);
    assert(neighbors[1].quad.level() == quad.level() + 1);
    assert(parent(neighbors[0].quad) == parent(neighbors[1].quad));
    const fub::p4est::quadrant<2> parent =
        fub::p4est::parent(neighbors[0].quad);
    assert(face_neighbor(quad, face) == parent);
    const int child_id_0 = child_id(neighbors[0].quad);
    const int child_id_1 = child_id(neighbors[1].quad);
    linearily_coarsen(grid.patch_data(neighbors[0].quad), interpolated,
                      child_id_0);
    linearily_coarsen(grid.patch_data(neighbors[1].quad), interpolated,
                      child_id_1);
    return interpolated;
  }

  // neighbor has a coarser level than quad
  assert(neighbors.size() == 1);
  if (neighbors[0].quad.level() == quad.level()) {
    return grid.patch_data(neighbors[0].quad);
  }

  assert(neighbors[0].quad.level() + 1 == quad.level());
  fub::p4est::quadrant<2> neighbor = face_neighbor(quad, face);
  assert(parent(neighbor) == neighbors[0].quad);
  const int child_id_0 = child_id(neighbor);
  linearily_refine(grid.patch_data(neighbors[0].quad), interpolated,
                   child_id_0);
  return interpolated;
}

void my_main(MPI_Comm communicator, int rank, int init_depth, int final_depth) {
  using namespace fub::v1::p4est;

  dynamic_extents_t<2> patch_extents(8, 8);
  uniform_cartesian_coordinates<2> coordinates({0.0, 0.0}, {1.0, 1.0},
                                               as_array(patch_extents));

  Grid::patch_data_info pdinfo{Variables(), patch_extents, coordinates};
  fub::p4est::connectivity<2> conn = unit_square;
  fub::p4est::forest<2> forest(communicator, conn, 0, init_depth, 1);
  fub::p4est::ghost_layer<2> ghost_layer(forest);
  Grid grid(std::move(forest), std::move(ghost_layer), pdinfo);

  span<const quadrant<2>> quadrants = grid.get_forest().trees()[0].quadrants();
  for (quadrant<2> quad : quadrants) {
    auto data = grid.patch_data(quad);
    auto coords = grid.coordinates(quad);
    my_initialize(data, coords);
  }
  fub::variable_data<Variables, double, 2> left_buffer(Variables(),
                                                       patch_extents);
  fub::variable_data<Variables, double, 2> right_buffer(Variables(),
                                                        patch_extents);

  for (int i = 0; i < final_depth - init_depth; ++i) {
    grid.exchange_ghost_data();

    fub::p4est::forest<2> forest = grid.get_forest();
    refine_if(forest, [&](int which_tree, quadrant<2> quad) {
      auto L2_derivatives_ = [&](fub::axis axis) {
        const face fl{axis, direction::left};
        const face fr{axis, direction::right};
        const std::array<double, 2> dx = grid.coordinates(quad).dx();
        Grid::const_patch_view middle = grid.patch_data(quad);
        Grid::const_patch_view left =
            interpolate_neighbor_data(grid, quad, fl, left_buffer.span());
        Grid::const_patch_view right =
            interpolate_neighbor_data(grid, quad, fr, right_buffer.span());
        const double volume = fub::accumulate(dx, 1.0, std::multiplies<>{});
        return L2_derivatives(as_int(axis), left, middle, right) > volume / 64;
      };
      return L2_derivatives_(axis::x) || L2_derivatives_(axis::y);
    });

    grid.interpolate_data_to(forest, [&](auto /* out */, auto in) {
      int size = in.quadrants.size();
      for (int i = 0; i < size; ++i) {
        my_initialize(in.datas[i], grid.coordinates(in.quadrants[i]));
      }
    });
  }

  struct State {
    Grid grid;
    std::chrono::duration<double> time;
    std::ptrdiff_t cycle;
  };
  using namespace std::literals;
  State state{std::move(grid), 1s, 1};
  output::cgns writer{};
  std::string filename = fmt::format("test.grid.p4est.{}.cgns", rank);
  auto file = writer.open(filename.c_str(), 2);
  writer.write(file, state);
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int init_depth = 0;
  if (argc > 1) {
    init_depth = std::atoi(argv[1]);
  }
  int final_depth = init_depth + 1;
  if (argc > 2) {
    final_depth = std::atoi(argv[2]);
  }
  my_main(MPI_COMM_WORLD, rank, init_depth, final_depth);
  MPI_Finalize();
}