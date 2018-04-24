// Copyright (c) 2017 Maikel Nadolski
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

#ifndef FUB_EULER_HLLESIGNALVELOCITY_HPP
#define FUB_EULER_HLLESIGNALVELOCITY_HPP

#include "fub/euler/ideal_gas.hpp"

namespace fub {
namespace euler {

struct hlle_signal_velocity {
public:
  using state_type = quantities<Density, Momentum<0>, SpeedOfSound>;

private:
  double m_beta{1.0};

public:
  explicit hlle_signal_velocity(double beta) : m_beta{beta} {};

  std::pair<double, double> do_computation(const state_type& left,
                                           const state_type& right) const
      noexcept;

  std::pair<scalar<double>, scalar<double>>
  do_simd_computation(simd_abi::scalar, const add_scalar_t<state_type>& left,
                      const add_scalar_t<state_type>& right) const noexcept;

#if defined(FUB_SIMD_HAS_SSE)
  std::pair<sse<double>, sse<double>> do_simd_computation(
      simd_abi::sse, const add_simd_t<state_type, simd_abi::sse>& left,
      const add_simd_t<state_type, simd_abi::sse>& right) const noexcept;
#endif

#if defined(FUB_SIMD_HAS_AVX)
  std::pair<avx<double>, avx<double>> do_simd_computation(
      simd_abi::avx, const add_simd_t<state_type, simd_abi::avx>& left,
      const add_simd_t<state_type, simd_abi::avx>& right) const noexcept;
#endif

public:
  template <typename M, int Rank>
  std::pair<double, double> operator()(const state_type& left,
                                       const state_type& right,
                                       const ideal_gas<M, Rank>&) const noexcept {
    return do_computation(left, right);
  }

  template <typename Abi, int Rank, typename L, typename R, typename M>
  auto operator()(const Abi& abi, const L& left, const R& right,
                  const ideal_gas<M, Rank>&) const noexcept {
    return do_simd_computation(abi, left, right);
  }
};

} // namespace euler
} // namespace fub

#endif // !HLLESIGNALVELOCITY_HPP
