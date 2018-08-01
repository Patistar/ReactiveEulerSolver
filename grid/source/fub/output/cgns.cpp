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

#include "fub/output/cgns.hpp"
#include "fub/layout_right.hpp"
#include "fub/tuple.hpp"

#include <array>
#include <numeric>
#include <vector>

namespace fub {
inline namespace v1 {
namespace output {
cgns_context::cgns_context(int file, int base) noexcept
    : m_file{file}, m_base{base} {}

cgns_context::cgns_context(cgns_context&& other) noexcept
    : m_file(std::exchange(other.m_file, -1)),
      m_base(std::exchange(other.m_base, -1)) {}

cgns_context& cgns_context::operator=(cgns_context&& other) noexcept {
  m_file = std::exchange(other.m_file, -1);
  m_base = std::exchange(other.m_base, -1);
  return *this;
}

cgns_context::~cgns_context() noexcept {
  if (m_file != -1) {
    cg_close(m_file);
  }
}

int cgns_context::file() const noexcept { return m_file; }
int cgns_context::base() const noexcept { return m_base; }

namespace {
const char* base_name = "Reactive Euler Solver";

[[noreturn]] void throw_cg_error() {
  std::string what = fmt::format(
      "cgns::ofile: An error occured while using the cgns library: {}.",
      cg_get_error());
  throw cgns_error(what);
}
} // namespace

namespace cgns_details {
void throw_if_cg_error(int return_value) {
  if (return_value) {
    throw_cg_error();
  }
}
} // namespace cgns_details

cgns_context cgns::open(const char* file_name, int dim) {
  int file;
  using cgns_details::throw_if_cg_error;
  throw_if_cg_error(cg_open(file_name, CG_MODE_WRITE, &file));
  int base;
  throw_if_cg_error(cg_base_write(file, base_name, dim, dim, &base));
  return cgns_context(file, base);
}

void cgns::iteration_data_write(const cgns_context& ctx,
                                std::chrono::duration<double> t,
                                long int cycle) {
  using cgns_details::throw_if_cg_error;
  throw_if_cg_error(
      cg_biter_write(ctx.file(), ctx.base(), "TimeIterValues", 1));
  throw_if_cg_error(
      cg_goto(ctx.file(), ctx.base(), "BaseIterativeData_t", 1, "end"));
  double t_ = t.count();
  cgsize_t size = 1;
  throw_if_cg_error(cg_array_write("TimeValues", RealDouble, 1, &size, &t_));
  throw_if_cg_error(
      cg_array_write("IterationValues", LongInteger, 1, &size, &cycle));
  throw_if_cg_error(
      cg_simulation_type_write(ctx.file(), ctx.base(), TimeAccurate));
}

namespace cgns_details {
namespace {

template <int Rank>
std::vector<double>
make_coordinates(const uniform_cartesian_coordinates<Rank>& coordinates,
                 int dim) {
  assert(0 <= dim && dim < Rank);
  auto extents = coordinates.extents();
  std::transform(extents.begin(), extents.end(), extents.begin(),
                 [](index e) { return e + 1; });
  const index size = std::accumulate(extents.begin(), extents.end(), index(1),
                                     std::multiplies<>());
  auto e = std::make_from_tuple<dynamic_extents_t<Rank>>(extents);
  layout_right::mapping<dynamic_extents_t<Rank>> mapping{e};
  std::vector<double> x(size);
  for_each_index(mapping, [&](auto... indices) {
    std::array<std::ptrdiff_t, Rank> reversed{indices...};
    std::reverse(reversed.begin(), reversed.end());
    std::array<double, Rank> xs = coordinates(indices...);
    x[fub::apply(mapping, reversed)] = xs[dim];
  });
  return x;
}

const char* coordinate_name(int dim) {
  switch (dim) {
  case 0:
    return "CoordinateX";
  case 1:
    return "CoordinateY";
  case 2:
    return "CoordinateZ";
  default:
    throw std::logic_error("Dimension is out of range.");
  }
}

template <std::size_t Rank>
zone_context zone_write_(const cgns_context& ctx, const char* name,
                         const std::array<index, Rank>& extents) {
  cgsize_t size[3][Rank];
  for (std::size_t dim = 0; dim < Rank; ++dim) {
    size[0][dim] = extents[dim] + 1; // Number of Vertices
    size[1][dim] = extents[dim];     // Number of Cells
    size[2][dim] = 0;                // This is always Zero for structured grids
  }
  int zone;
  throw_if_cg_error(cg_zone_write(ctx.file(), ctx.base(), name, &size[0][0],
                                  Structured, &zone));
  return {ctx.file(), ctx.base(), zone};
}

template <>
zone_context zone_write_<std::size_t(1)>(const cgns_context& ctx,
                                         const char* name,
                                         const std::array<index, 1>& extents) {
  cgsize_t size[3][2];
  size[0][0] = extents[0] + 1; // Number of Vertices
  size[1][0] = extents[0];     // Number of Cells
  size[2][0] = 0;              // This is always Zero for structured grids
  size[0][1] = 2;              // Number of Vertices
  size[1][1] = 1;              // Number of Cells
  size[2][1] = 0;              // This is always Zero for structured grids
  int zone;
  throw_if_cg_error(cg_zone_write(ctx.file(), ctx.base(), name, &size[0][0],
                                  Structured, &zone));
  return {ctx.file(), ctx.base(), zone};
}

/// @brief Writes coordinates wrt a specified octant into a specified cgns file.
template <int Rank>
void coordinates_write_(const zone_context& zone,
                        const uniform_cartesian_coordinates<Rank>& coords) {
  for (int dim = 0; dim < Rank; ++dim) {
    // Create and write the coordinate values
    {
      std::vector<double> x = make_coordinates<Rank>(coords, dim);
      int index_coord = 0;
      throw_if_cg_error(cg_coord_write(zone.file, zone.base, zone.id,
                                       RealDouble, coordinate_name(dim),
                                       x.data(), &index_coord));
    }
    // Add Dimensional Info for this coordinate direction. (Its given in meters)
    throw_if_cg_error(cg_goto(zone.file, zone.base, "Zone_t", zone.id,
                              "GridCoordinates", 0, coordinate_name(dim), 0,
                              "end"));
    throw_if_cg_error(cg_dataclass_write(Dimensional));
    throw_if_cg_error(cg_units_write(Kilogram, Meter, Second, Kelvin, Radian));
    double exponents[]{0, 1, 0, 0, 0};
    throw_if_cg_error(cg_exponents_write(RealDouble, &exponents[0]));
  }
}

template <>
void coordinates_write_<1>(const zone_context& zone,
                           const uniform_cartesian_coordinates<1>& coords) {
  std::array<double, 2> lower{coords.lower()[0], 0};
  std::array<double, 2> upper{coords.upper()[0], 1};
  std::array<index, 2> extents{coords.extents()[0], 1};
  uniform_cartesian_coordinates<2> lifted{lower, upper, extents};
  coordinates_write_(zone, lifted);
}
} // namespace

void coordinates_write(const zone_context& zone,
                       const uniform_cartesian_coordinates<1>& coordinates) {
  return coordinates_write_(zone, coordinates);
}
void coordinates_write(const zone_context& zone,
                       const uniform_cartesian_coordinates<2>& coordinates) {
  return coordinates_write_(zone, coordinates);
}
void coordinates_write(const zone_context& zone,
                       const uniform_cartesian_coordinates<3>& coordinates) {
  return coordinates_write_(zone, coordinates);
}

zone_context zone_write(const cgns_context& ctx, const char* name,
                        const std::array<index, 1>& extents) {
  return zone_write_(ctx, name, extents);
}
zone_context zone_write(const cgns_context& ctx, const char* name,
                        const std::array<index, 2>& extents) {
  return zone_write_(ctx, name, extents);
}
zone_context zone_write(const cgns_context& ctx, const char* name,
                        const std::array<index, 3>& extents) {
  return zone_write_(ctx, name, extents);
}

int cell_centered_field_write(const zone_context& zone, int solution,
                              const char* name, span<const double> view) {
  int field_number;
  throw_if_cg_error(cg_field_write(zone.file, zone.base, zone.id, solution,
                                   RealDouble, name, view.data(),
                                   &field_number));
  return field_number;
}

} // namespace cgns_details
} // namespace output
} // namespace v1
} // namespace fub
