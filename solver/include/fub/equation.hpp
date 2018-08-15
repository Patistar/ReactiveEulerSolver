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

#ifndef FUB_SOLVER_EQUATION_HPP
#define FUB_SOLVER_EQUATION_HPP

#include "fub/mdspan.hpp"
#include "fub/variable.hpp"
#include "fub/quantities.hpp"
#include "fub/variable_list.hpp"
#include "fub/variable_view.hpp"

namespace fub {
inline namespace v1 {
// Precision Type

template <typename Eq> using floating_point_t = typename Eq::floating_point;

// Span Types
template <typename Eq> using span_t = span<floating_point_t<Eq>>;

template <typename Eq>
using mdspan_t =
    basic_mdspan<floating_point_t<Eq>, dynamic_extents_t<Eq::rank()>>;

template <typename Eq>
using const_mdspan_t =
    basic_mdspan<floating_point_t<Eq>, dynamic_extents_t<Eq::rank()>>;

// Variable Lists

template <typename Eq> using complete_t = typename Eq::complete;

template <typename Eq> using conservative_t = typename Eq::conservative;

// Patche Types

template <typename Eq>
using state_view_t = variable_view<complete_t<Eq>, floating_point_t<Eq>,
                                   dynamic_extents_t<Eq::rank()>>;

template <typename Eq>
using patch_t = variable_view<complete_t<Eq>, floating_point_t<Eq>,
                              dynamic_extents_t<Eq::rank()>>;

template <typename Eq>
using const_patch_t = variable_view<complete_t<Eq>, const floating_point_t<Eq>,
                                    dynamic_extents_t<Eq::rank()>>;

template <typename Eq>
using fluxes_t = variable_view<conservative_t<Eq>, floating_point_t<Eq>,
                               dynamic_extents_t<Eq::rank()>>;

template <typename Eq>
using const_fluxes_t =
    variable_view<conservative_t<Eq>, const floating_point_t<Eq>,
                  dynamic_extents_t<Eq::rank()>>;

// Reference Types

template <typename Eq, typename A = accessor_native<floating_point_t<Eq>>>
using state_ref_t = variable_ref<complete_t<Eq>, floating_point_t<Eq>, A>;

template <typename Eq, typename A = accessor_native<const floating_point_t<Eq>>>
using const_state_ref_t =
    variable_ref<complete_t<Eq>, const floating_point_t<Eq>, A>;

template <typename Eq, typename A = accessor_native<floating_point_t<Eq>>>
using flux_ref_t = variable_ref<conservative_t<Eq>, floating_point_t<Eq>, A>;

template <typename Eq, typename A = accessor_native<const floating_point_t<Eq>>>
using const_flux_ref_t =
    variable_ref<conservative_t<Eq>, const floating_point_t<Eq>, A>;

template <typename Eq, typename A = accessor_native<floating_point_t<Eq>>>
using conservatives_ref_t =
    variable_ref<conservative_t<Eq>, floating_point_t<Eq>, A>;

template <typename Eq, typename A = accessor_native<const floating_point_t<Eq>>>
using const_conservatives_ref_t =
    variable_ref<conservative_t<Eq>, const floating_point_t<Eq>, A>;

// State Types

template <typename Eq>
using flux_state_t =
    basic_variable_data<conservative_t<Eq>,
                        mdspan<floating_point_t<Eq>, 1>>;


} // namespace v1
} // namespace fub

#endif