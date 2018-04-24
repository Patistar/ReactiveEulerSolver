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

#ifndef FUB_FLUXMETHOD_HPP
#define FUB_FLUXMETHOD_HPP

#include "fub/simd.hpp"
#include "fub/type_traits.hpp"
#include "fub/variables.hpp"
#include <type_traits>

namespace fub {

template <typename S, typename L, typename R, typename E>
using has_compute_numeric_flux_t =
    decltype(std::declval<const S&>().compute_numeric_flux(
        std::declval<const L&>(), std::declval<const R&>(),
        std::declval<const E&>()));

template <typename RiemannSolver, typename Equation>
struct is_simd_enabled
    : is_detected<has_compute_numeric_flux_t, RiemannSolver,
                  add_simd_t<typename Equation::CompleteState>,
                  add_simd_t<typename Equation::CompleteState>, Equation> {};

} // namespace fub

#endif // !FLUXMETHOD_HPP
