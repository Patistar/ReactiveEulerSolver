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

#ifndef FUB_EQUATION_HPP
#define FUB_EQUATION_HPP

#include "fub/variables.hpp"

namespace fub {

template <typename Eq> struct complete_state {
  using type = typename Eq::complete_state;
};
template <typename Eq>
using complete_state_t = typename complete_state<Eq>::type;

template <typename Eq> struct conservative_state {
  using type = typename Eq::conservative_state;
};
template <typename Eq>
using conservative_state_t = typename conservative_state<Eq>::type;

template <typename Eq>
using conservative_variables_t = as_tuple_t<conservative_state_t<Eq>>;

template <typename Eq> struct flux_type {
  using type = add_flux_t<conservative_state_t<Eq>>;
};
template <typename Eq> using flux_type_t = typename flux_type<Eq>::type;

template <typename Eq> struct species_tuple_ {
  using type = typename Eq::species_tuple;
};
template <typename Eq> using species_tuple_t = typename species_tuple_<Eq>::type;

template <typename Eq> species_tuple_t<Eq> species_tuple(const Eq&) noexcept {
  return {};
}

} // namespace fub

#endif
