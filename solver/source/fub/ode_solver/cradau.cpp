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

#include <blaze/Math.h>

#include "fub/ode_solver/cradau.hpp"
#include "fub/algorithm.hpp"
#include "fub/optional.hpp"


namespace fub {
namespace ode_solver {
namespace {

using ode_system_t = radau5::ode_system_t;
using feedback_t = radau5::feedback_t;

enum class radau5_code { accepted, rejected, restart, error };

using matrix_type = blaze::CustomMatrix<double, blaze::unaligned,
                                        blaze::unpadded, blaze::columnMajor>;

using vector_type =
    blaze::CustomVector<double, blaze::unaligned, blaze::unpadded>;

struct integration_state {
  span<char> available_memory;
  matrix_type jacobian;
  matrix_type e1;
  matrix_type e2;
  vector_type w;
  vector_type z;
  vector_type a;
  vector_type f_y;
  vector_type scaling;
  span<int> ipiv1;
  span<int> ipiv2;
  double dw_norm_old{0.};
  double time{0.};
  double step_size{0.};
  double step_size_old{1.};
  double step_size_inv{1.};
  double theta_old{1.};
  double tolerance{1.};
  double error{1.};
  double error_old{1.};
  double stopping_criterium;
  index step_count{0};
  int restart_count{0};
  int newton_iteration_count{0};
};

/// \brief Returns the time coefficients for the 3-stage RADAU IIa solver.
///
/// These are most often denoted as the `c` vector in the `(A, b, c)`-butcher
/// scheme.
constexpr std::array<double, 3> get_time_coefficients() noexcept {
  constexpr double sq6 =
      2.449489742783178098197284074705891391965947480656670128432;
  return {{(4. - sq6) / 10., (4. + sq6) / 10., 1.}};
}

blaze::StaticMatrix<double, 3, 3> get_transform_matrix() noexcept {
  static constexpr double T11 = 9.1232394870892942792E-02;
  static constexpr double T12 = -0.14125529502095420843E0;
  static constexpr double T13 = -3.0029194105147424492E-02;
  static constexpr double T21 = 0.24171793270710701896E0;
  static constexpr double T22 = 0.20412935229379993199E0;
  static constexpr double T23 = 0.38294211275726193779E0;
  static constexpr double T31 = 0.96604818261509293619E0;
  static constexpr double T32 = 1.0;
  static constexpr double T33 = 0.0;
  return blaze::StaticMatrix<double, 3, 3>{
      {T11, T12, T13}, {T21, T22, T23}, {T31, T32, T33}};
}

blaze::StaticMatrix<double, 3, 3> get_transform_inverse_matrix() noexcept {
  static constexpr double TI11 = 4.3255798900631553510E0;
  static constexpr double TI12 = 0.33919925181580986954E0;
  static constexpr double TI13 = 0.54177053993587487119E0;
  static constexpr double TI21 = -4.1787185915519047273E0;
  static constexpr double TI22 = -0.32768282076106238708E0;
  static constexpr double TI23 = 0.47662355450055045196E0;
  static constexpr double TI31 = -0.50287263494578687595E0;
  static constexpr double TI32 = 2.5719269498556054292E0;
  static constexpr double TI33 = -0.59603920482822492497E0;
  return blaze::StaticMatrix<double, 3, 3>{
      {TI11, TI12, TI13}, {TI21, TI22, TI23}, {TI31, TI32, TI33}};
}

blaze::StaticMatrix<double, 3, 3> get_lambda_matrix() noexcept {
  // std::pow(81., 1. / 3.)
  static constexpr double tr81 =
      4.326748710922225146964914932340328765175607760498051732639;
  // std::pow(9., 1. / 3.)
  static constexpr double tr9 =
      2.080083823051904114530056824357885386337805340373262109697;
  // std::sqrt(3)
  static constexpr double sq3 =
      1.732050807568877293527446341505872366942805253810380628055;

  // Compute 'Eigen' values
  static constexpr double eigenvalue = 30. / (6. + tr81 - tr9);
  static constexpr double alpha_ = (12. - tr81 + tr9) / 60.;
  static constexpr double beta_ = (tr81 + tr9) * sq3 / 60.;
  static constexpr double norm2 = alpha_ * alpha_ + beta_ * beta_;
  static constexpr double alpha = alpha_ / norm2;
  static constexpr double beta = beta_ / norm2;

  // Matrix
  static constexpr double L11 = eigenvalue;
  static constexpr double L12 = 0.0;
  static constexpr double L13 = 0.0;
  static constexpr double L21 = 0.0;
  static constexpr double L22 = alpha;
  static constexpr double L23 = -beta;
  static constexpr double L31 = 0.0;
  static constexpr double L32 = beta;
  static constexpr double L33 = alpha;

  return blaze::StaticMatrix<double, 3, 3>{
      {L11, L12, L13}, {L21, L22, L23}, {L31, L32, L33}};
}

span<char> drop_until_aligned(span<char> memory, index alignment) {
  const intptr_t address = reinterpret_cast<intptr_t>(memory.data());
  const intptr_t offset = address % alignment;
  const intptr_t shift = alignment - offset;
  return drop(memory, shift);
}

span<char> consume_memory(integration_state& state, index size,
                          index alignment) {
  span<char> aligned_memory =
      drop_until_aligned(state.available_memory, alignment);
  span<char> memory = take(aligned_memory, size);
  state.available_memory = drop(aligned_memory, size);
  return memory;
}

matrix_type new_matrix(integration_state& state, index rows, index cols) {
  index size = rows * cols;
  index byte_size = sizeof(double) * size;
  span<char> memory = consume_memory(state, byte_size, alignof(double));
  double* data = new (memory.data()) double[size];
  return matrix_type(data, rows, cols);
}

vector_type new_vector(integration_state& state, index size) {
  index byte_size = sizeof(double) * size;
  span<char> memory = consume_memory(state, byte_size, alignof(double));
  double* data = new (memory.data()) double[size];
  return vector_type(data, size);
}

span<int> new_ipiv(integration_state& state, index size) {
  index byte_size = sizeof(int) * size;
  span<char> memory = consume_memory(state, byte_size, alignof(int));
  int* data = new (memory.data()) int[size];
  return span<int>(data, size);
}

void allocate_subobjects(integration_state& state, index size) {
  constexpr int stages = 3;
  state.jacobian = new_matrix(state, size, size);
  state.e1 = new_matrix(state, size, size);
  state.e2 = new_matrix(state, 2 * size, 2 * size);
  state.a = new_vector(state, stages * size);
  state.w = new_vector(state, stages * size);
  state.z = new_vector(state, stages * size);
  state.scaling = new_vector(state, size);
  state.f_y = new_vector(state, size);
  state.ipiv1 = new_ipiv(state, size);
  state.ipiv2 = new_ipiv(state, size);
}

template <typename Left, typename Right, typename Out>
void kronecker_product(const Left& left, const Right& right, Out& out) {
  for (std::size_t i = 0; i < blaze::rows(left); ++i) {
    for (std::size_t j = 0; j < blaze::columns(left); ++j) {
      const std::size_t row_lower = i * blaze::rows(right);
      const std::size_t col_lower = j * blaze::columns(right);
      auto block = blaze::submatrix(out, row_lower, col_lower,
                                    blaze::rows(right), blaze::columns(right));
      block = left(i, j) * right;
    }
  }
}

void assemble_iteration_matrix(integration_state& state,
                               const radau5& /* options */) {
  const blaze::StaticMatrix<double, 3, 3> lambda = get_lambda_matrix();
  const std::size_t n_rows = state.jacobian.rows();
  const std::size_t n_cols = state.jacobian.columns();
  assert(n_rows == n_cols);
  const blaze::IdentityMatrix<double> identity(n_rows);
  // Assemble E1
  const double h_ = state.step_size_inv;
  state.e1 = h_ * lambda(0, 0) * identity - state.jacobian;
  blaze::getrf(state.e1, state.ipiv1.data());
  // Assemble E2
  kronecker_product(h_ * submatrix(lambda, 1, 1, 2, 2), identity, state.e2);
  submatrix(state.e2, 0, 0, n_rows, n_cols) -= state.jacobian;
  submatrix(state.e2, n_rows, n_cols, n_rows, n_cols) -= state.jacobian;
  blaze::getrf(state.e2, state.ipiv2.data());
}

void solve_iteration_matrix(integration_state& state,
                            const radau5& /* options */,
                            ode_system_t ode_system, const vector_type& y_0) {
  static constexpr auto time_coefficients = get_time_coefficients();
  auto lambda = get_lambda_matrix();
  const auto T_inv = get_transform_inverse_matrix();
  index size = state.jacobian.rows();

  auto a_0 = subvector(state.a, 0, size);
  auto a_1 = subvector(state.a, size, size);
  auto a_2 = subvector(state.a, 2 * size, size);
  auto z_0 = subvector(state.z, 0, size);
  auto z_1 = subvector(state.z, size, size);
  auto z_2 = subvector(state.z, 2 * size, size);

  // Compute F(Z^k) and store it into the `z`-vector

  a_0 = y_0 + z_0;
  a_1 = y_0 + z_1;
  a_2 = y_0 + z_2;
  const double dt = state.step_size;
  ode_system(state.time + time_coefficients[0] * dt, a_0, z_0);
  ode_system(state.time + time_coefficients[1] * dt, a_1, z_1);
  ode_system(state.time + time_coefficients[2] * dt, a_2, z_2);

  // Compute (T_inv o Id) F(Z^k) and store it into the `z`-vector

  a_0 = z_0;
  a_1 = z_1;
  a_2 = z_2;
  z_0 = T_inv(0, 0) * a_0 + T_inv(0, 1) * a_1 + T_inv(0, 2) * a_2;
  z_1 = T_inv(1, 0) * a_0 + T_inv(1, 1) * a_1 + T_inv(1, 2) * a_2;
  z_2 = T_inv(2, 0) * a_0 + T_inv(2, 1) * a_1 + T_inv(2, 2) * a_2;

  // Compute `(h_inv Lambda o Id) W^k - (T_inv o Id) F(Z^k)` and store it in `a`

  const auto w_0 = subvector(state.w, 0, size);
  const auto w_1 = subvector(state.w, size, size);
  const auto w_2 = subvector(state.w, 2 * size, size);
  double h_inv = state.step_size_inv;
  a_0 = z_0 - h_inv * lambda(0, 0) * w_0;
  a_1 = z_1 - h_inv * (lambda(1, 1) * w_1 + lambda(1, 2) * w_2);
  a_2 = z_2 - h_inv * (lambda(2, 1) * w_1 + lambda(2, 2) * w_2);

  // Solve the linear system and compute DW^k, store it in `a`.

  auto a_12 = subvector(state.a, size, 2 * size);

  blaze::getrs(state.e1, a_0, 'N', state.ipiv1.data());
  blaze::getrs(state.e2, a_12, 'N', state.ipiv2.data());
}

void update_z_and_w(integration_state& state) {
  state.w += state.a;

  const auto T = get_transform_matrix();
  index size = state.jacobian.rows();

  const auto w_0 = subvector(state.w, 0, size);
  const auto w_1 = subvector(state.w, size, size);
  const auto w_2 = subvector(state.w, 2 * size, size);

  auto z_0 = subvector(state.z, 0, size);
  auto z_1 = subvector(state.z, size, size);
  auto z_2 = subvector(state.z, 2 * size, size);

  z_0 = T(0, 0) * w_0 + T(0, 1) * w_1 + T(0, 2) * w_2;
  z_1 = T(1, 0) * w_0 + T(1, 1) * w_1 + T(1, 2) * w_2;
  z_2 = T(2, 0) * w_0 + T(2, 1) * w_1 + T(2, 2) * w_2;
}

radau5_code do_first_newton_iteration_step(integration_state& state,
                                           const radau5& options,
                                           ode_system_t ode_system,
                                           const vector_type& y) {
  solve_iteration_matrix(state, options, ode_system, y);
  const index size = y.size();
  auto a_0 = subvector(state.a, 0, size);
  auto a_1 = subvector(state.a, size, size);
  auto a_2 = subvector(state.a, 2 * size, size);
  state.dw_norm_old =
      std::sqrt((sqrNorm(a_0 / state.scaling) + sqrNorm(a_1 / state.scaling) +
                 sqrNorm(a_2 / state.scaling)) /
                3 * size);
  state.theta_old = 1;
  update_z_and_w(state);
  state.newton_iteration_count += 1;
  if (state.dw_norm_old < 1e-12) {
    return radau5_code::accepted;
  } else {
    return radau5_code::rejected;
  }
}

radau5_code do_newton_iteration_step(integration_state& state,
                                     const radau5& options,
                                     ode_system_t ode_system,
                                     const vector_type& y_0) noexcept {
  solve_iteration_matrix(state, options, ode_system, y_0);
  const index size = y_0.size();
  auto a_0 = subvector(state.a, 0, size);
  auto a_1 = subvector(state.a, size, size);
  auto a_2 = subvector(state.a, 2 * size, size);
  double dw_norm_new =
      std::sqrt((sqrNorm(a_0 / state.scaling) + sqrNorm(a_1 / state.scaling) +
                 sqrNorm(a_2 / state.scaling)) /
                3 * size);
  double theta = dw_norm_new / state.dw_norm_old;
  state.dw_norm_old = dw_norm_new;
  state.theta_old = std::exchange(theta, std::sqrt(theta * state.theta_old));
  if (theta < 1.0) {
    const double eta = theta / (1 - theta);
    const int k_max = options.max_newton_iteration_count;
    const int k = state.newton_iteration_count;
    const double theta_max = std::pow(theta, k_max - k - 1) / k;
    const double expected_iteration_error = eta * theta_max * dw_norm_new;
    if (expected_iteration_error < state.stopping_criterium) {
      state.newton_iteration_count += 1;
      update_z_and_w(state);
      const double iteration_error = eta * dw_norm_new;
      return (iteration_error < state.stopping_criterium)
                 ? radau5_code::accepted
                 : radau5_code::rejected;
    }
  }
  state.step_size *= 0.5;
  state.step_size_inv = 1.0 / state.step_size;
  assemble_iteration_matrix(state, options);
  state.restart_count += 1;
  state.newton_iteration_count = 0;
  state.theta_old = 1;
  blaze::reset(state.w);
  blaze::reset(state.z);
  return radau5_code::restart;
}

radau5_code do_newton_iteration(integration_state& state, const radau5& options,
                                ode_system_t ode_system,
                                const vector_type& y_0) {
  radau5_code code{};
  do {
    code = do_newton_iteration_step(state, options, ode_system, y_0);
  } while (code == radau5_code::rejected);
  return code;
}

void compute_jacobian(integration_state& state, const radau5& options,
                      ode_system_t ode_system, const vector_type& y_0) {
  constexpr double eps = std::numeric_limits<double>::epsilon();
  const index size = y_0.size();
  auto dy = subvector(state.a, 0, size);
  auto dydx = subvector(state.a, size, size);
  ode_system(state.time, y_0, state.f_y);
  dy = y_0;
  for (index i = 0; i < size; ++i) {
    const double delta = std::sqrt(eps * std::max(1.0E-5, std::fabs(y_0[i])));
    dy[i] += delta;
    ode_system(state.time, dy, dydx);
    column(state.jacobian, i) = (dydx - state.f_y) / delta;
    dy[i] = y_0[i];
  }
  assemble_iteration_matrix(state, options);
}

void estimate_error(integration_state& state, const radau5&,
                    ode_system_t ode_system, const vector_type& y) {
  // std::sqrt(6)
  static constexpr double sq6 =
      2.449489742783178098197284074705891391965947480656670128432;

  // See [Hairer] p.123
  static constexpr std::array<double, 3> e{
      {(-13. - 7. * sq6) / 3., (-13. + 7. * sq6) / 3., -1. / 3.}};

  const index size = y.size();

  const double h = state.step_size;

  const auto z_0 = subvector(state.z, 0, size);
  const auto z_1 = subvector(state.z, size, size);
  const auto z_2 = subvector(state.z, 2 * size, size);

  auto error = subvector(state.a, 0, size);
  auto term = subvector(state.a, size, size);

  term = (e[0] * z_0 + e[1] * z_1 + e[2] * z_2) / h;
  error = state.f_y + term;
  blaze::getrs(state.e1, error, 'N', state.ipiv1.data());
  double error_norm =
      std::min(std::sqrt(sqrNorm(error / state.scaling) / size), 1.0E-10);
  if (error_norm > 1.0) {
    error += y;
    auto f_error = subvector(state.a, 2 * size, size);
    ode_system(state.time, error, f_error);
    blaze::getrs(state.e1, f_error, 'N', state.ipiv1.data());
    f_error += term;
    error_norm =
        std::min(std::sqrt(sqrNorm(error / state.scaling) / size), 1.0E-10);
  }
  state.error_old = std::exchange(state.error, error_norm);
}

void update_step_size(integration_state& state, const radau5& options) {
  const index k_max = options.max_newton_iteration_count;
  const index k = state.newton_iteration_count;
  const double fac = 0.9 * double(2 * k_max - 1) / double(2 * k - 1);
  const double err4 = std::sqrt(std::sqrt(state.error));
  const double h = state.step_size;
  const double h_new = [&] {
    if (state.step_count > 1) {
      const double err4_old = std::sqrt(std::sqrt(state.error_old));
      const double h_old = state.step_size_old;
      return fac * (h / err4) * (h / h_old) * (err4_old / err4);
    } else {
      return fac * h / err4;
    }
  }();
  if (h_new < h || 1.2 * h < h_new) {
    // Enforce: 0.2 <= (h_new / h) <= 8.0
    const double h_limited = clamp(h_new, 0.2 * h, 8. * h);
    state.step_size_old = std::exchange(state.step_size, h_limited);
    state.step_size_inv = 1.0 / state.step_size;
    assemble_iteration_matrix(state, options);
  }
}

void update_scaling(integration_state& state, const radau5& options,
                  const vector_type& y) {
  for (std::size_t i = 0; i < y.size(); ++i) {
    state.scaling[i] = options.absolute_tolerance +
                       options.relative_tolerance * std::abs(y[i]);
  }
}

void advance_in_time(integration_state& state, const radau5& options,
                     ode_system_t ode_system, vector_type& y) {
  estimate_error(state, options, ode_system, y);
  index size = y.size();
  y += subvector(state.z, 2 * size, size);
  state.time += state.step_size;
  state.step_count += 1;
  update_step_size(state, options);
  state.restart_count = 0;
  update_scaling(state, options, y);
}

radau5_code try_step(integration_state& state, const radau5& options,
                     ode_system_t ode_system, vector_type& y) {
  radau5_code code{};
  do {
    code = do_first_newton_iteration_step(state, options, ode_system, y);
    if (code != radau5_code::accepted) {
      code = do_newton_iteration(state, options, ode_system, y);
    }
  } while (code == radau5_code::restart);
  assert(code == radau5_code::accepted || code == radau5_code::error);
  if (code == radau5_code::accepted) {
    advance_in_time(state, options, ode_system, y);
  }
  return code;
}

void invoke_feedback(const optional<feedback_t>& feedback,
                     const integration_state& state, const vector_type& y) {
  if (feedback) {
    radau5::integration_info info{};
    info.newton_iteration_count = state.newton_iteration_count;
    (*feedback)(state.time, span<const double>(y.data(), y.size()), info);
  }
}

void integrate_impl(integration_state& state, const radau5& options,
                    optional<feedback_t> feedback, ode_system_t ode_system,
                    vector_type& y, double dt) {
  constexpr double eps = std::numeric_limits<double>::epsilon();
  radau5_code code{};
  while (state.time < dt && code != radau5_code::error) {
    invoke_feedback(feedback, state, y);
    const double upper_bound = dt - state.time + 10 * eps;
    state.step_size = std::min(state.step_size, upper_bound);
    if (state.newton_iteration_count != 1 && 1.0E-3 < state.theta_old) {
      compute_jacobian(state, options, ode_system, y);
    }
    code = try_step(state, options, ode_system, y);
  }
  invoke_feedback(feedback, state, y);
}

} // namespace

void radau5::integrate(ode_system_t ode_system, std::array<double, 2> x,
                       span<double> y_0, span<char> memory,
                       optional<feedback_t> feedback) const {
  integration_state state;
  state.available_memory = memory;
  state.step_size = initial_step_size;
  state.step_size_inv = 1.0 / state.step_size;
  state.theta_old = 1.0;
  state.stopping_criterium = std::min(0.03, std::sqrt(relative_tolerance));
  vector_type y(y_0.data(), y_0.size());
  allocate_subobjects(state, y.size());
  update_scaling(state, *this, y);
  integrate_impl(state, *this, std::move(feedback), ode_system, y, x[1] - x[0]);
}

} // namespace ode_solver
} // namespace fub