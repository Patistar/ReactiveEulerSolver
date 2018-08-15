#include "fub/euler/perfect_gas.hpp"
#include "fub/utility.hpp"
#include "fub/variable_data.hpp"

#include <benchmark/benchmark.h>

#include <random>

constexpr int Rank = 3;

using Equation = fub::euler::perfect_gas<Rank>;
using Equationf = fub::euler::perfect_gas<Rank, float>;
using Complete = Equation::complete;
using Conservative = Equation::conservative;
using Extents = fub::dynamic_extents_t<Rank>;

template <typename S> using View = Equation::view<S>;
template <typename S> using ConstView = Equation::const_view<S>;

template <typename S> using Viewf = Equationf::view<S>;
template <typename S> using ConstViewf = Equationf::const_view<S>;

template <typename S, typename A>
using Patch = fub::variable_data<S, A, Rank>;

void reconstruct_double(Equation equation, ConstView<Conservative> cons,
                        View<Complete> patch);

void reconstruct_double_simd(Equation equation, ConstView<Conservative> cons,
                             View<Complete> patch);

void reconstruct_float(Equationf equation, ConstViewf<Conservative> cons,
                       Viewf<Complete> patch);

void reconstruct_float_simd(Equationf equation, ConstViewf<Conservative> cons,
                            Viewf<Complete> patch);

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
    q[momentum<1, 0>] = rho * u;
    q[energy] = rho * (e_int + e_kin);
  }
}

template <int Rank> Extents make_extents(std::ptrdiff_t size) {
  std::array<std::ptrdiff_t, Rank> array;
  array.fill(size);
  return fub::apply([](auto... extent) { return Extents(extent...); }, array);
}

static void BM_reconstruct_double(benchmark::State& state) {
  Patch<Complete, double> patch(Complete(), make_extents<Rank>(state.range(0)));
  Patch<Conservative, double> cons(Conservative(),
                                   make_extents<Rank>(state.range(0)));
  Equation equation{1.4, 28.};
  random_initialize(equation, cons);
  for (auto _ : state) {
    reconstruct_double(equation, fub::as_const(cons), patch);
  }
}
BENCHMARK(BM_reconstruct_double)->RangeMultiplier(2)->Range(4, 4 << 4);

static void BM_reconstruct_float(benchmark::State& state) {
  Patch<Complete, float> patch(Complete(), make_extents<Rank>(state.range(0)));
  Patch<Conservative, float> cons(Conservative(),
                                  make_extents<Rank>(state.range(0)));
  Equationf equation{1.4, 28.};
  random_initialize(equation, cons);
  for (auto _ : state) {
    reconstruct_float(equation, fub::as_const(cons), patch);
  }
}
BENCHMARK(BM_reconstruct_float)->RangeMultiplier(2)->Range(4, 4 << 4);

static void BM_reconstruct_double_simd(benchmark::State& state) {
  Patch<Complete, double> patch(Complete(), make_extents<Rank>(state.range(0)));
  Patch<Conservative, double> cons(Conservative(),
                                   make_extents<Rank>(state.range(0)));
  Equation equation{1.4, 28.};
  random_initialize(equation, cons);
  for (auto _ : state) {
    reconstruct_double_simd(equation, fub::as_const(cons), patch);
  }
}
BENCHMARK(BM_reconstruct_double_simd)->RangeMultiplier(2)->Range(4, 4 << 4);

static void BM_reconstruct_float_simd(benchmark::State& state) {
  Patch<Complete, float> patch(Complete(), make_extents<Rank>(state.range(0)));
  Patch<Conservative, float> cons(Conservative(),
                                  make_extents<Rank>(state.range(0)));
  Equationf equation{1.4, 28.};
  random_initialize(equation, cons);
  for (auto _ : state) {
    reconstruct_float_simd(equation, fub::as_const(cons), patch);
  }
}
BENCHMARK(BM_reconstruct_float_simd)->RangeMultiplier(2)->Range(4, 4 << 4);

BENCHMARK_MAIN();