#include "fub/euler/perfect_gas.hpp"

constexpr int Rank = 3;

using Equation = fub::euler::perfect_gas<Rank>;
using Complete = Equation::complete;
using Conservative = Equation::conservative;
using Extents = fub::dynamic_extents_t<Rank>;

template <typename S> using View = Equation::view<S>;
template <typename S> using ConstView = Equation::const_view<S>;

using Equationf = fub::euler::perfect_gas<Rank, float>;
using Completef = Equationf::complete;
using Conservativef = Equationf::conservative;
template <typename S> using Viewf = Equationf::view<S>;
template <typename S> using ConstViewf = Equationf::const_view<S>;

void flux_double(Equation equation, ConstView<Complete> patch,
                 View<Conservative> fluxes) {
  equation.flux(patch, fluxes);
}

void flux_float(Equationf equation, ConstViewf<Completef> patch,
                Viewf<Conservativef> fluxes) {
  equation.flux(patch, fluxes);
}

void reconstruct_double(Equation equation, ConstView<Conservative> cons,
                        View<Complete> patch) {
  equation.reconstruct(cons, patch);
}

void reconstruct_float(Equationf equation, ConstViewf<Conservativef> cons,
                       Viewf<Completef> patch) {
  equation.reconstruct(cons, patch);
}

void reconstruct_double_simd(Equation equation, ConstView<Conservative> cons,
                             View<Complete> patch) {
  equation.reconstruct(fub::execution::vec, cons, patch);
}

void reconstruct_float_simd(Equationf equation, ConstViewf<Conservativef> cons,
                             Viewf<Completef> patch) {
  equation.reconstruct(fub::execution::vec, cons, patch);
}