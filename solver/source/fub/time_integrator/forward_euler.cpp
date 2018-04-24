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

#include "fub/time_integrator/forward_euler.hpp"

namespace fub {
namespace time_integrator {
namespace {
template <typename T>
auto do_integrate_impl_(T lambda, T u, T fL, T fR) noexcept {
  return u + lambda * (fL - fR);
}
} // namespace
double forward_euler::do_integrate_impl(double lambda, double u, double fL,
                                       double fR) noexcept {
  return do_integrate_impl_(lambda, u, fL, fR);
}

scalar<double> forward_euler::do_integrate_impl(simd_abi::scalar,
                                               scalar<double> lambda,
                                               scalar<double> u,
                                               scalar<double> fL,
                                               scalar<double> fR) noexcept {
  return do_integrate_impl_(lambda, u, fL, fR);
}

#if defined(FUB_SIMD_HAS_SSE)
simd<double> forward_euler::do_integrate_impl(simd_abi::native<double>,
                                             simd<double> lambda,
                                             simd<double> u, simd<double> fL,
                                             simd<double> fR) noexcept {
  return do_integrate_impl_(lambda, u, fL, fR);
}
#endif
} // namespace time_integrator
} // namespace fub
