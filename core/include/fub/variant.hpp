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

/// @file This header selects an avialable implementation of std::variant and
/// makes it available throgh the `fub::` namespace.

#ifndef FUB_CORE_VARIANT_HPP
#define FUB_CORE_VARIANT_HPP

#include "fub/core/config.hpp"

#if defined(FUB_WITH_STD_VARIANT)
#include <variant>
namespace fub {
using std::get;
using std::get_if;
using std::in_place_index;
using std::variant;
using std::visit;
} // namespace fub
#elif __has_include(<mpark/variant.hpp>)
#include <mpark/variant.hpp>
namespace fub {
using mpark::get;
using mpark::get_if;
using mpark::in_place_index;
using mpark::variant;
using mpark::visit;
} // namespace fub
#else
#error("No implementation for std::variant could be found.")
#endif

#endif // !FUB_CORE_VARIANT_HPP
