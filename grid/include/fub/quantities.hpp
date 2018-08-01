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

#ifndef FUB_GRID_QUANTITIES_HPP
#define FUB_GRID_QUANTITIES_HPP

#include "fub/variable_data.hpp"

namespace fub {
inline namespace v1 {

template <typename VariableList, typename T>
class quantities : private basic_variable_data<VariableList, T, mdspan<T, 1>> {
private:
  using Base_ = basic_variable_data<VariableList, T, mdspan<T, 1>>;

public:
  quantities() = default;

  template <typename... Args,
            typename = std::enable_if_t<sizeof...(Args) ==
                                        VariableList::static_size()>,
            typename = std::enable_if_t<
                conjunction<std::is_convertible<Args, T>...>::value>>
  quantities(Args... data) : Base_() {
    auto datas = hana::make_tuple(data...);
    hana::for_each(
        hana::make_range(hana::size_c<0>, hana::size_c<sizeof...(Args)>),
        [&](auto Is) { span()[Is] = datas[Is]; });
  }

  template <typename Tag> T& operator[](Tag tag) {
    return Base_::operator[](tag)(0);
  }

  template <typename Tag> const T& operator[](Tag tag) const {
    return Base_::operator[](tag)(0);
  }

  auto span() noexcept { return Base_::span(); }
  auto span() const noexcept { return Base_::span(); }
};

} // namespace v1
} // namespace fub

#endif