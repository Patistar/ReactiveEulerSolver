#include "fub/euler/hlle_signal_velocity.hpp"
#include "fub/euler/perfect_gas.hpp"
#include "fub/utility.hpp"
#include "fub/variable_data.hpp"

#include <benchmark/benchmark.h>

#include <random>

constexpr int Rank = 3;

using Equation = fub::euler::perfect_gas<Rank>;
using EquationF = fub::euler::perfect_gas<Rank, float>;
using EquationLD = fub::euler::perfect_gas<Rank, long double>;
using Complete = Equation::complete;
using Conservative = Equation::conservative;
using Extents = fub::dynamic_extents_t<Rank>;

template <typename S> using View = Equation::view<S>;
template <typename S> using ConstView = Equation::const_view<S>;

template <typename S> using ViewF = EquationF::view<S>;
template <typename S> using ConstViewF = EquationF::const_view<S>;

template <typename S, typename A> using Patch = fub::variable_data<S, A, Rank>;

fub::euler::hlle_signal_velocities<double>
signal_double(Equation equation, fub::const_state_ref_t<Equation> qL,
              fub::const_state_ref_t<Equation> qR);

fub::euler::hlle_signal_velocities<float>
signal_float(EquationF equation, fub::const_state_ref_t<EquationF> qL,
             fub::const_state_ref_t<EquationF> qR);

template <typename Eq>
void random_initialize(Eq equation,
                       typename Eq::template view<Conservative> cons) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> rho_dis(0.1, 5.0);
  std::uniform_real_distribution<> u_dis(-100., 100.);
  std::uniform_real_distribution<> T_dis(100., 1000.);
  const auto cv = equation.heat_capacity_at_constant_volume();
  for (auto q : cons) {
    using namespace fub::euler;
    const double rho = rho_dis(gen);
    const double u = u_dis(gen);
    const double T = T_dis(gen);
    const double e_int = cv * T;
    const double e_kin = 0.5 * u * u * rho;
    q[density] = rho;
    q[momentum<Eq::rank()>] = rho * u;
    q[energy] = rho * (e_int + e_kin);
  }
}

template <int Rank> Extents make_extents(std::ptrdiff_t size) {
  std::array<std::ptrdiff_t, Rank> array;
  array.fill(size);
  return fub::apply([](auto... extent) { return Extents(extent...); }, array);
}

static void BM_hlle_signal_double(benchmark::State& state) {
  Patch<Complete, double> patch(Complete(), make_extents<Rank>(state.range(0)));
  Patch<Conservative, double> cons(Conservative(),
                                   make_extents<Rank>(state.range(0)));
  std::vector<double> lefts(patch.size());
  std::vector<double> rights(patch.size());
  fub::mdspan_t<Equation> left_view(lefts, patch.get_extents());
  fub::mdspan_t<Equation> right_view(rights, patch.get_extents());
  Equation equation{1.4, 28.};
  random_initialize(equation, cons);
  fub::for_each(
      [=](const fub::const_conservatives_ref_t<Equation>& q,
          const fub::state_ref_t<Equation>& w) { equation.reconstruct(q, w); },
      cons, patch);
  for (auto _ : state) {
    auto rows = fub::rows(patch);
    auto left_signal = left_view.span().begin();
    auto right_signal = right_view.span().begin();
    for (auto row : rows) {
      auto left = row.begin();
      auto right = std::next(left);
      auto end = row.end();
      while (right != end) {
        auto signals = signal_double(equation, *left, *right);
        *left_signal++ = signals.left;
        *right_signal++ = signals.right;
        ++left;
        ++right;
      }
    }
  }
}
BENCHMARK(BM_hlle_signal_double)->RangeMultiplier(2)->Range(4, 4 << 4);

static void BM_hlle_signal_float(benchmark::State& state) {
  Patch<Complete, float> patch(Complete(), make_extents<Rank>(state.range(0)));
  Patch<Conservative, float> cons(Conservative(),
                                  make_extents<Rank>(state.range(0)));
  std::vector<float> lefts(patch.size());
  std::vector<float> rights(patch.size());
  fub::mdspan_t<EquationF> left_view(lefts, patch.get_extents());
  fub::mdspan_t<EquationF> right_view(rights, patch.get_extents());
  EquationF equation{1.4, 28.};
  random_initialize(equation, cons);
  fub::for_each(
      [=](const fub::const_conservatives_ref_t<EquationF>& q,
          const fub::state_ref_t<EquationF>& w) { equation.reconstruct(q, w); },
      cons, patch);
  for (auto _ : state) {
    auto rows = fub::rows(patch);
    auto left_signal = left_view.span().begin();
    auto right_signal = right_view.span().begin();
    for (auto row : rows) {
      auto left = row.begin();
      auto right = std::next(left);
      auto end = row.end();
      while (right != end) {
        auto signals = signal_float(equation, *left, *right);
        *left_signal++ = signals.left;
        *right_signal++ = signals.right;
        ++left;
        ++right;
      }
    }
  }
}
BENCHMARK(BM_hlle_signal_float)->RangeMultiplier(2)->Range(4, 4 << 4);

BENCHMARK_MAIN();