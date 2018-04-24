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

#include "fub/euler/hlle_signal_velocity.hpp"

#include <algorithm>
#include <cmath>
#include <type_traits>
#include <utility>

namespace fub {
namespace euler {

std::pair<double, double>
hlle_signal_velocity::do_computation(const state_type& left,
                                     const state_type& right) const noexcept {
  using namespace variables;
  using std::sqrt;
  const auto sqRhoL = sqrt(left[density]);
  const auto sqRhoR = sqrt(right[density]);
  const auto uL = left[momentum<0>] / left[density];
  const auto uR = right[momentum<0>] / right[density];
  const auto aL = left[speed_of_sound];
  const auto aR = right[speed_of_sound];
  const auto roeU = (sqRhoL * uL + sqRhoR * uR) / (sqRhoL + sqRhoR);
  const auto roeA = (sqRhoL * aL + sqRhoR * aR) / (sqRhoL + sqRhoR) +
                    0.5f * (sqRhoL * sqRhoR) / (sqRhoL + sqRhoR) *
                        (sqRhoL + sqRhoR) * (uR - uL) * (uR - uL);
  const auto sL1 = uL - m_beta * aL;
  const auto sL2 = roeU - roeA;
  const auto sR1 = roeU + roeA;
  const auto sR2 = uR + m_beta * aR;
  return std::make_pair(std::min(sL1, sL2), std::max(sR1, sR2));
}

namespace {
template <typename Abi>
std::pair<simd<double, Abi>, simd<double, Abi>>
impl(const Abi&, double beta,
     const add_simd_t<hlle_signal_velocity::state_type, Abi>& left,
     const add_simd_t<hlle_signal_velocity::state_type, Abi>& right) noexcept {
  using namespace variables;
  const auto sqRhoL = sqrt(left[density]);
  const auto sqRhoR = sqrt(right[density]);
  const auto uL = left[momentum<0>] / left[density];
  const auto uR = right[momentum<0>] / right[density];
  const auto aL = left[speed_of_sound];
  const auto aR = right[speed_of_sound];
  const auto roeU = (sqRhoL * uL + sqRhoR * uR) / (sqRhoL + sqRhoR);
  const auto roeA = (sqRhoL * aL + sqRhoR * aR) / (sqRhoL + sqRhoR) +
                    0.5f * (sqRhoL * sqRhoR) / (sqRhoL + sqRhoR) *
                        (sqRhoL + sqRhoR) * (uR - uL) * (uR - uL);
  const auto sL1 = uL - beta * aL;
  const auto sL2 = roeU - roeA;
  const auto sR1 = roeU + roeA;
  const auto sR2 = uR + beta * aR;
  return std::make_pair(min(sL1, sL2), max(sR1, sR2));
}
} // namespace

std::pair<scalar<double>, scalar<double>>
hlle_signal_velocity::do_simd_computation(
    simd_abi::scalar, const add_scalar_t<state_type>& left,
    const add_scalar_t<state_type>& right) const noexcept {
  return impl(simd_abi::scalar(), m_beta, left, right);
}

#if defined(FUB_SIMD_HAS_SSE)
std::pair<sse<double>, sse<double>> hlle_signal_velocity::do_simd_computation(
    simd_abi::sse, const add_simd_t<state_type, simd_abi::sse>& left,
    const add_simd_t<state_type, simd_abi::sse>& right) const noexcept {
  return impl(simd_abi::sse(), m_beta, left, right);
}
#endif

#if defined(FUB_SIMD_HAS_AVX)
std::pair<avx<double>, avx<double>> hlle_signal_velocity::do_simd_computation(
    simd_abi::avx, const add_simd_t<state_type, simd_abi::avx>& left,
    const add_simd_t<state_type, simd_abi::avx>& right) const noexcept {
  return impl(simd_abi::avx(), m_beta, left, right);
}
#endif

} // namespace euler
} // namespace fub
