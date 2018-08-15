#include "fub/euler/hlle_riemann_solver.hpp"
#include "fub/euler/perfect_gas.hpp"
#include "fub/godunov_method.hpp"
#include "fub/time_integrator/forward_euler.hpp"
#include "fub/variable_data.hpp"

#include <benchmark/benchmark.h>

#include <random>

constexpr int Rank = 3;

using Equation = fub::euler::perfect_gas<Rank>;
using Extents = fub::dynamic_extents_t<Rank>;
using FluxMethod = fub::godunov_method<fub::euler::hlle_riemann_solver>;
template <typename S, typename A> using Patch = fub::variable_data<S, A, Rank>;

void integrate_double(const Equation& equation, const FluxMethod& method,
                      std::chrono::duration<double> dt, double dx,
                      fub::const_patch_t<Equation> left,
                      fub::const_patch_t<Equation> mid,
                      fub::const_patch_t<Equation> right,
                      fub::patch_t<Equation> next,
                      fub::fluxes_t<Equation> fluxes);

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

static void BM_forward_euler_integrate_double(benchmark::State& state) {
  Extents extents = make_extents<Rank>(state.range(0));
  FluxMethod godunov{};
  std::chrono::duration<double> dt(0.1);
  double dx = 0.1;
  Patch<fub::complete_t<Equation>, double> left(fub::complete_t<Equation>(),
                                                extents);
  
  Patch<fub::complete_t<Equation>, double> middle(fub::complete_t<Equation>(),
                                                  extents);
  
  Patch<fub::complete_t<Equation>, double> right(fub::complete_t<Equation>(),
                                                 extents);
  
  Patch<fub::complete_t<Equation>, double> next(fub::complete_t<Equation>(),
                                                extents);
  
  Patch<fub::conservative_t<Equation>, double> initial(
      fub::conservative_t<Equation>(), extents);
  
  Patch<fub::conservative_t<Equation>, double> fluxes(
      fub::conservative_t<Equation>(), fub::grow<0>(extents));

  Equation equation{1.4, 28.};

  random_initialize(equation, initial, left);

  random_initialize(equation, initial, middle);

  random_initialize(equation, initial, right);
  
  for (auto _ : state) {
    integrate_double(equation, godunov, dt, dx, left, middle, right, next,
                     fluxes);
  }
}
BENCHMARK(BM_forward_euler_integrate_double)
    ->RangeMultiplier(2)
    ->Range(4, 4 << 4);
BENCHMARK_MAIN();