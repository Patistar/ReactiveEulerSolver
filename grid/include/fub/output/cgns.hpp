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

#ifndef FUB_OUTPUT_CGNS_HPP
#define FUB_OUTPUT_CGNS_HPP

#include "fub/uniform_cartesian_coordinates.hpp"
#include "fub/variable_view.hpp"

extern "C" {
#include <cgnslib.h>
}

#include <chrono>
#include <string>

#include <fmt/format.h>

namespace fub {
inline namespace v1 {
namespace output {
/// This is the base error class for exceptions thrown by this module.
///
/// Every exception which is thrown by any cgns related function is inherited
/// from this base class.
struct cgns_error : public std::runtime_error {
  using std::runtime_error::runtime_error;
};

/// This wraps a cgns file descriptor and base id number which are used in cgns
/// functions.
///
/// This class owns its file descriptor. Thus the file gets closed if this
/// object gets destroyed.
class cgns_context;

struct cgns {
  /// Opens a CGNS file and returns a cgns_context pointing to that file.
  ///
  /// This function also creates a Base in the CGNS file with the specified
  /// problem dimension.
  ///
  /// \return Returns a cgns_context object which refers to an open file and
  /// cgns base.
  static cgns_context open(const char* file_name, int dimension, int mode = CG_MODE_WRITE);

  template <typename Grid, typename Box>
  static void write(const cgns_context& ctx, const Grid& grid, const Box& box);

  template <typename State>
  static void write(const cgns_context& ctx, const State& state);

  static void iteration_data_write(const cgns_context& ctx,
                                   std::chrono::duration<double> t,
                                   long int cycle);
};

/// This wraps a cgns file descriptor and base id number which are used in cgns
/// functions.
///
/// This class owns its file descriptor. Thus the file gets closed if this
/// object gets destroyed.
class cgns_context {
public:
  cgns_context() = delete;
  cgns_context(const cgns_context&) = delete;
  cgns_context& operator=(const cgns_context&) = delete;
  cgns_context(cgns_context&& other) noexcept;
  cgns_context& operator=(cgns_context&& other) noexcept;

  /// Closes the file referred by the file descriptor.
  ~cgns_context() noexcept;

  /// Returns the file descriptor.
  int file() const noexcept;

  /// Returns the cgns base id number.
  int base() const noexcept;

private:
  friend struct cgns;

  cgns_context(int file, int base) noexcept;

  int m_file;
  int m_base;
};

namespace cgns_details {
void throw_if_cg_error(int return_value);

struct zone_context {
  int file;
  int base;
  int id;
};

/// @{
zone_context zone_write(const cgns_context& context, const char* name,
                        const std::array<index, 1>& extents);

zone_context zone_write(const cgns_context& context, const char* name,
                        const std::array<index, 2>& extents);

zone_context zone_write(const cgns_context& context, const char* name,
                        const std::array<index, 3>& extents);
/// @}

/// @brief Fills a specified CGNS zone with structured coordinates.
/// @{
void coordinates_write(const zone_context& zone,
                       const uniform_cartesian_coordinates<1>&);

void coordinates_write(const zone_context& zone,
                       const uniform_cartesian_coordinates<2>&);

void coordinates_write(const zone_context& zone,
                       const uniform_cartesian_coordinates<3>&);
/// @}

/// @brief Writes a data std::array into the CGNS file.
int cell_centered_field_write(const zone_context& zone, int solution,
                              const char* name, span<const double>);
/// @}

template <typename Patch>
void flow_solution_write(const zone_context& zone, const Patch& patch) {
  int solution;
  throw_if_cg_error(cg_sol_write(zone.file, zone.base, zone.id, "FlowSolution",
                                 CellCenter, &solution));
  for_each(patch.get_variable_list(), [&](auto var) {
    cell_centered_field_write(zone, solution, var.name(), patch[var].span());
  });
}

template <typename Quad> std::array<char, 33> make_zone_name(Quad quad) {
  std::array<char, 33> array{};
  fmt::basic_memory_buffer<char, 33> buffer;
  auto x = quad.coordinates();
  fmt::format_to(buffer, "Zone-{}_{}_{}", quad.level(), x[0], x[1]);
  std::copy_n(buffer.data(), buffer.size(), array.data());
  return array;
}
}; // namespace cgns_details

template <typename Grid, typename Box>
void cgns::write(const cgns_context& ctx, const Grid& grid, const Box& box) {
  using namespace cgns_details;
  std::array<char, 33> zone_name = cgns_details::make_zone_name(box);
  auto coordinates = grid.coordinates(box);
  auto patch = grid.patch_data(box);
  zone_context zone = zone_write(ctx, zone_name.data(), coordinates.extents());
  coordinates_write(zone, coordinates);
  flow_solution_write(zone, patch);
}

template <typename State>
void cgns::write(const cgns_context& ctx, const State& state) {
  iteration_data_write(ctx, state.time, state.cycle);
  auto trees = state.grid.get_forest().trees();
  for (auto&& tree : trees) {
    auto quads = tree.quadrants();
    for (auto quad : quads) {
      write(ctx, state.grid, quad);
    }
  }
}

} // namespace output
} // namespace v1
} // namespace fub

#endif // !CGNS_HPP
