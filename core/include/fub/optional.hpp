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

/// @file This header selects an avialable implementation of std::optional and
/// makes it available throgh the `fub::` namespace.

#ifndef FUB_CORE_OPTIONAL_HPP
#define FUB_CORE_OPTIONAL_HPP

#include "fub/core/config.hpp"

#if defined(FUB_WITH_STD_OPTIONAL)
#include <optional>
#elif defined(FUB_WITH_STD_EXPERIMENTAL_OPTIONAL)
#include <experimental/optional>
#elif __has_include("nonstd/optional.hpp")
#include "nonstd/optional.hpp"
// #elif __has_include(<boost/optional.hpp>)
// #include <boost/optional.hpp>
#else
#error("No implementation for std::optional could be found.")
#endif

namespace fub {

#if defined(FUB_WITH_STD_OPTIONAL)
using std::nullopt;
using std::optional;
#elif defined(FUB_WITH_STD_EXPERIMENTAL_OPTIONAL)
using std::experimental::nullopt;
using std::experimental::optional;
#elif __has_include("nonstd/optional.hpp")
// #else
using nonstd::optional;
using nonstd::nullopt;
// #elif __has_include(<boost/optional.hpp>)
// using boost::optional;
// static const auto nullopt = boost::none;
#endif

} // namespace fub

#endif
