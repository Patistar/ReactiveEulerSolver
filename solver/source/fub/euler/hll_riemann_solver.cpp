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

#include "fub/euler/hll_riemann_solver.hpp"
#include <cassert>

namespace fub {
namespace euler {
/// @brief Returns the original HLL flux for a single specified quantity.
///
/// @todo citation
double hll_riemann_solver::apply_formula(double uL, double uR, double fL,
                                         double fR, double sL,
                                         double sR) noexcept {
  assert(sL < sR);
  if (0 <= sL) {
    return fL;
  } else if (sR <= 0) {
    return fR;
  }
  return (sR * fL - sL * fR + sL * sR * (uR - uL)) / (sR - sL);
}

namespace {
template <typename Abi>
simd<double, Abi>
apply_formula_impl(const simd<double, Abi>& uL, const simd<double, Abi>& uR,
                   const simd<double, Abi>& fL, const simd<double, Abi>& fR,
                   const simd<double, Abi>& sL,
                   const simd<double, Abi>& sR) noexcept {
  simd<double, Abi> result =
      (sR * fL - sL * fR + sL * sR * (uR - uL)) / (sR - sL);
  where(sL >= 0, result) = fL;
  where(sL < 0 && sR <= 0, result) = fR;
  return result;
}
} // namespace

scalar<double> hll_riemann_solver::apply_formula(
    simd_abi::scalar, const scalar<double>& uL, const scalar<double>& uR,
    const scalar<double>& fL, const scalar<double>& fR,
    const scalar<double>& sL, const scalar<double>& sR) noexcept {
  return apply_formula_impl(uL, uR, fL, fR, sL, sR);
}

#if defined(FUB_SIMD_HAS_SSE)
sse<double>
hll_riemann_solver::apply_formula(simd_abi::sse, const sse<double>& uL,
                                  const sse<double>& uR, const sse<double>& fL,
                                  const sse<double>& fR, const sse<double>& sL,
                                  const sse<double>& sR) noexcept {
  return apply_formula_impl(uL, uR, fL, fR, sL, sR);
}
#endif

#if defined(FUB_SIMD_HAS_AVX)
avx<double>
hll_riemann_solver::apply_formula(simd_abi::avx, const avx<double>& uL,
                                  const avx<double>& uR, const avx<double>& fL,
                                  const avx<double>& fR, const avx<double>& sL,
                                  const avx<double>& sR) noexcept {
  return apply_formula_impl(uL, uR, fL, fR, sL, sR);
}
#endif

} // namespace euler
} // namespace fub
