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

#ifndef FUB_EULER_HLLESIGNALVELOCITY_HPP
#define FUB_EULER_HLLESIGNALVELOCITY_HPP

#include "fub/equation.hpp"
#include "fub/euler/variables.hpp"
#include "fub/simd.hpp"

namespace fub {
inline namespace v1 {
namespace euler {

template <typename T> struct hlle_signal_velocities {
  T left;
  T right;
};

/// A signal velocity policy which can be used with approximate Riemann Solvers
/// of the HLL type.
class hlle_signal_velocity {
public:
  template <typename Equation, typename StateRef>
  hlle_signal_velocities<typename StateRef::value_type>
  operator()(Equation equation, const StateRef& qL, const StateRef& qR) const {
    using float_v = typename StateRef::value_type;
    using std::max;
    using std::min;
    using std::sqrt;
    constexpr int Rank = Equation::rank();
    const float_v rhoL = equation.get(density, qL);
    const float_v rhoR = equation.get(density, qR);
    assert(all_of(rhoL > 0.0 || rhoR > 0.0));
    const float_v sqRhoL = sqrt(rhoL);
    const float_v sqRhoR = sqrt(rhoR);
    const float_v uL = equation.get(velocity<Rank>, qL);
    const float_v uR = equation.get(velocity<Rank>, qR);
    const float_v aL = equation.get(speed_of_sound, qL);
    const float_v aR = equation.get(speed_of_sound, qR);
    assert(all_of(aL >= 0.0));
    assert(all_of(aR >= 0.0));
    const float_v roeU = (sqRhoL * uL + sqRhoR * uR) / (sqRhoL + sqRhoR);
    const float_v roeA =
        sqrt((sqRhoL * aL * aL + sqRhoR * aR * aR) / (sqRhoL + sqRhoR) +
             (sqRhoL * sqRhoR) / 2 / ((sqRhoL + sqRhoR) * (sqRhoL + sqRhoR)) *
                 (uR - uL) * (uR - uL));
    const float_v sL1 = uL - aL;
    const float_v sL2 = roeU - roeA;
    const float_v sR1 = roeU + roeA;
    const float_v sR2 = uR + aR;
    assert(all_of(min(sL1, sL2) <= max(sR1, sR2)));
    return {min(sL1, sL2), max(sR1, sR2)};
  }
};

} // namespace euler
} // namespace v1
} // namespace fub

#endif // !HLLESIGNALVELOCITY_HPP
