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

#ifndef FUB_SOLVER_GRID_HPP
#define FUB_SOLVER_GRID_HPP

#include "fub/type_traits.hpp"
#include "fub/variable_data.hpp"

namespace fub {
inline namespace v1 {

template <typename Grid> using quadrant_t = typename remove_cvref_t<Grid>::quadrant_type;

template <typename Grid>
using patch_data_t =
    decltype(std::declval<Grid>().patch_data(quadrant_t<Grid>{}));

template <typename Grid>
using coordinates_t = decltype(std::declval<Grid>().coordinates());

template <typename Grid>
using variable_list_t = decltype(std::declval<Grid>().variable_list());

template <typename Grid, typename VariableList = variable_list_t<Grid>>
using patch_buffer_t = variable_data<VariableList, typename Grid::value_type, Grid::rank()>; 

} // namespace v1
} // namespace fub

#endif