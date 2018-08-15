#include "fub/equation.hpp"
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

void solve_double(Equation equation, fub::const_patch_t<Equation> left,
                  fub::const_patch_t<Equation> middle,
                  fub::const_patch_t<Equation> right,
                  fub::fluxes_t<Equation> flux);

template <typename Eq>
void random_initialize(Eq equation, fub::fluxes_t<Eq> cons,
                       fub::patch_t<Eq> patch) {
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
  fub::for_each(
      [=](const fub::const_conservatives_ref_t<Equation>& q,
          const fub::state_ref_t<Equation>& w) { equation.reconstruct(q, w); },
      cons, patch);
}

template <int Rank> Extents make_extents(std::ptrdiff_t size) {
  std::array<std::ptrdiff_t, Rank> array;
  array.fill(size);
  return fub::apply([](auto... extent) { return Extents(extent...); }, array);
}

static void BM_hlle_solve_double(benchmark::State& state) {
  Extents extents = make_extents<Rank>(state.range(0));
  Patch<Complete, double> left(Complete(), extents);
  Patch<Complete, double> middle(Complete(), extents);
  Patch<Complete, double> right(Complete(), extents);
  Patch<Conservative, double> initial(Conservative(), extents);
  Patch<Conservative, double> fluxes(Conservative(), fub::grow<0>(extents));
  Equation equation{1.4, 28.};
  random_initialize(equation, initial, left);
  random_initialize(equation, initial, middle);
  random_initialize(equation, initial, right);
  for (auto _ : state) {
    solve_double(equation, left, middle, right, fluxes);
  }
}
BENCHMARK(BM_hlle_solve_double)->RangeMultiplier(2)->Range(4, 4 << 4);

BENCHMARK_MAIN();