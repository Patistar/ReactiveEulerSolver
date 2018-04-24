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

#include "fub/interval_map.hpp"
#include <iostream>
#include <iomanip>
#include <string>

template <typename Map>
void print(std::ostream& out, const Map& map) {
  for (auto&& [key, mapped] : map) {
    out << '{' << key << ", " << mapped << "}\n";
  }
  out << '\n';
}

int main() {
  fub::interval_map<std::uint64_t, std::uint64_t> map{};

  print(std::cout, map.get_map());
  map.insert(1, 3, 1);
  map.insert(3, 5, 2);
  map.insert(4, 7, 3);
  print(std::cout, map.get_map());

  map.insert(0, 7, 1);
  map.insert(0, 3, 2);
  print(std::cout, map.get_map());
}
