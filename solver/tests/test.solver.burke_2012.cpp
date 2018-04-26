// Copyright (c) 2018 Maikel Nadolski
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "fub/euler/mechanism/burke_2012.hpp"
#include "fub/euler/kinetic_source_term.hpp"

#include <iostream>
#include <vector>

int main()
{
  using namespace fub::euler::mechanism::burke2012;
  const int n_species = std::tuple_size<Burke2012::species_tuple>::value;
  std::vector<double> x(n_species);
  const std::ptrdiff_t h2 = Index_v<H2>;
  const std::ptrdiff_t o2 = Index_v<O2>;
  x[Index_v<H2>] = 2.0;
  x[Index_v<O2>] = 1.0;
  fub::euler::ideal_gas<Burke2012> eq{};
  auto state = eq.set_TPX(1100, 4e5, x);
  auto source_term = fub::euler::make_kinetic_source_term(eq);
  using namespace std::chrono_literals;
  auto next = source_term.advance_state(state, 1.e-3s);
  fub::for_each_tuple_element([&](auto var) {
    std::cout << var.name() << ": " << next[var] << '\n';
  }, get_variables(next));
}
