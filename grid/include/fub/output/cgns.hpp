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

#include "fub/array.hpp"
#include "fub/octree.hpp"
#include "fub/patch_view.hpp"
#include "fub/span.hpp"
#include "fub/uniform_cartesian_coordinates.hpp"

extern "C" {
#include <cgnslib.h>
}

#include <chrono>
#include <string>

namespace fub {
namespace output {
struct cgns_error : public std::runtime_error {
  using std::runtime_error::runtime_error;
};

class cgns_context {
public:
  cgns_context() = delete;
  cgns_context(const cgns_context&) = delete;
  cgns_context& operator=(const cgns_context&) = delete;
  cgns_context(cgns_context&& other) noexcept;
  cgns_context& operator=(cgns_context&& other) noexcept;
  cgns_context(int file, int base) noexcept;
  ~cgns_context() noexcept;

  int file() const noexcept;
  int base() const noexcept;

private:
  int m_file;
  int m_base;
};

struct cgns {
  /// @brief Opens a CGNS file and returns a cgns_context pointing to that file.
  ///
  /// This function also creates a Base in the CGNS file with the specified
  /// problem dimension.
  static cgns_context open(const char* file_name, int dimension);

  template <int Rank, typename View, typename Equation>
  static void write(const cgns_context& ctx, const octant<Rank>& octant,
                    const View& patch,
                    const uniform_cartesian_coordinates<Rank>& coordinates,
                    const Equation& equation);

  template <int Rank, typename Grid, typename Equation>
  static void write_grid(cgns_context& ctx, const Grid& grid,
                         const uniform_cartesian_coordinates<Rank>& coordinates,
                         const Equation& equation);

  static void iteration_data_write(const cgns_context& ctx,
                                   std::chrono::duration<double> t,
                                   long int cycle);
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
                        const array<index, 1>& extents);

zone_context zone_write(const cgns_context& context, const char* name,
                        const array<index, 2>& extents);

zone_context zone_write(const cgns_context& context, const char* name,
                        const array<index, 3>& extents);
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

/// @brief Writes a data array into the CGNS file.
int cell_centered_field_write(const zone_context& zone, int solution,
                              const char* name, span<const double>);

/// @brief This fills the specified zone with the user-provided dimensional
/// information.
void physical_dimension_write_impl(const zone_context& zone, int solution,
                                   const char* variable,
                                   const array<double, 5>&);

/// @brief We mark a specified variable as NonDimensional data.
void physical_dimension_write_impl(const zone_context& zone, int solution,
                                   const char* variable);

/// @brief Here `var.physical_dimensions()` is a valid expression and we use it
/// to get the dimensional exponents.
template <typename Variable>
void physical_dimension_write_(const zone_context& zone, int solution,
                               Variable var, std::true_type) {
  physical_dimension_write_impl(zone, solution, var.name(),
                                var.physical_dimensions());
}

/// @brief Here `var.physical_dimensions()` is *not* a valid expression and we
/// mark the variable as non-dimensional.
template <typename Variable>
void physical_dimension_write_(const zone_context& zone, int solution,
                               Variable var, std::false_type) {
  physical_dimension_write_impl(zone, solution, var.name());
}

/// @brief This function fills dimensional information for the specified
/// variable.
template <typename Variable>
void physical_dimension_write(const zone_context& zone, int solution,
                              Variable var) {
  physical_dimension_write_(zone, solution, var,
                            variable_is_dimensional<Variable>());
}
/// @}

template <typename View, typename Equation, typename... Vars>
void flow_solution_write(const zone_context& zone, const View& patch,
                         const Equation& equation,
                         const std::tuple<Vars...>& vars) {
  int solution;
  throw_if_cg_error(cg_sol_write(zone.file, zone.base, zone.id, "FlowSolution",
                                 CellCenter, &solution));
  std::vector<double> buffer(patch.extents().size());
  span<const double> view{buffer.data(), static_cast<index>(buffer.size())};
  for_each_tuple_element(
      [&](auto var) {
        span<double> out{buffer.data(), static_cast<index>(buffer.size())};
        for_each_row(
            [&](auto row) {
              std::transform(
                  row.begin(), row.end(), out.begin(),
                  [&](const auto& q) { return equation.get(var, q); });
              out = span<double>{out.data() + row.size(),
                                 out.size() - row.size()};
            },
            patch);
        cell_centered_field_write(zone, solution, var.name(), view);
        // physical_dimension_write(zone, solution, var);
      },
      vars);
}

array<char, 33> make_zone_name(const octant<1>& o) noexcept;
array<char, 33> make_zone_name(const octant<2>& o) noexcept;
array<char, 33> make_zone_name(const octant<3>& o) noexcept;

}; // namespace cgns_details

template <int Rank, typename View, typename Equation>
void cgns::write(const cgns_context& ctx, const octant<Rank>& octant,
                 const View& patch,
                 const uniform_cartesian_coordinates<Rank>& coordinates,
                 const Equation& equation) {
  using namespace cgns_details;
  array<char, 33> zone_name = make_zone_name(octant);
  zone_context zone = zone_write(ctx, zone_name.data(), coordinates.extents());
  coordinates_write(zone, adapt(coordinates, octant));
  flow_solution_write(zone, patch, equation, view_variables_t<View>());
}

} // namespace output
} // namespace fub

#endif // !CGNS_HPP
