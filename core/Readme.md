# Core library

## Introduction

This is a small C++14 compatible collection of library types and functions which 
are only available for more recent C++ standards like C++17 or C++20. These 
include

### Class Templates

 * `fub::extents<Es...>` - A compact multi-dimensional integral cartesian product.
 * `fub::optional<T>` - Depending on availability it is either `std::optional<T>`, `std::experimental::optional<T>` or `boost::optional<T>`.
 * `fub::variant<T1, ..., Tn>` - Depending on availability it is either `std::variant<T1, ... Tn>` or `mpark::variant<T1, ..., Tn>`.
 * `fub::span<T, N>` - A one-dimensinoal view over a contiguous array of values.  
 * `fub::mdspan<T, extents<Es...>>` - A multi-dimensional view over a contiguous array.
 * `fub::simd<T, Abi>` - A templated type wrapper for intrinsic vector operations. 
 
### Function Templates
 * `fub::apply` - A replacement for `std::apply`.
 * `fub::invoke` - A replacement for `std::invoke`.
 * `fub::visit` - A replacement for `std::visit`.
 * `fub::clamp` - A replacement for `std::clamp`.
 * `fub::accumulate` - A constexpr-enabled version `std::accumulate`.
 * `fub::transform_reduce` - A constexpr-enabled verison `std::transform_reduce`.
 * `fub::for_each_tuple_element` - Invoke a function object for each element of a tuple.
 * `fub::foldl` - Left fold over all tuple elements.

### Helper Traits
 * `fub::is_detected`
 * `fub::conjunction`
 * `fub::disjunction`
 * `fub::negation`
 
### Type List Manipulation

 * `fub::head`
 * `fub::tail`
 * `fub::concat`