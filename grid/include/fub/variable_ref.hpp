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

#ifndef FUB_GRID_VARIABLE_REF_HPP
#define FUB_GRID_VARIABLE_REF_HPP

#include "fub/variable_list.hpp"
#include <array>

namespace fub {
inline namespace v1 {
template <typename VariableView> class variable_ref {
public:
  static constexpr index rank = VariableView::extents_type::rank();

  variable_ref() = delete;
  variable_ref(const variable_ref&) = delete;

  variable_ref(VariableView& view,
               const std::array<index, rank>& indices) noexcept
      : m_view{view}, m_indices{indices} {}

  template <typename Variables>
  variable_ref& operator=(const Variables& variables) {
    for_each(get_variable_list(), [&](auto variable) {
      fub::apply(m_view[variable], m_indices) = variables[variable];
    });
    return *this;
  }

  template <typename Tag>
  constexpr decltype(auto) operator[](Tag tag) noexcept {
    return fub::apply(m_view[tag], m_indices);
  }

  template <typename Tag>
  constexpr decltype(auto) operator[](Tag tag) const noexcept {
    return fub::apply(m_view[tag], m_indices);
  }

  constexpr decltype(auto) get_variable_list() const noexcept {
    return m_view.get_variable_list();
  }

private:
  VariableView& m_view;
  std::array<index, rank> m_indices;
};
} // namespace v1
} // namespace fub

#endif