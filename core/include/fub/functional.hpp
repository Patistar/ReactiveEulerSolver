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

#ifndef FUB_CORE_FUNCTIONAL_HPP
#define FUB_CORE_FUNCTIONAL_HPP

#include "fub/type_traits.hpp"

#ifdef FUB_WITH_STD_INVOKE
#include <functional>
namespace fub {
inline namespace v1 { using std::invoke; }
} // namespace fub
#else
#include <range/v3/utility/invoke.hpp>
namespace fub {
inline namespace v1 { using ranges::invoke; }
} // namespace fub
#endif

namespace fub {
inline namespace v1 {

template <typename> class function_ref;

template <typename R, typename... Args> class function_ref<R(Args...)> {
public:
  function_ref(R (*function_pointer)(Args...))
      : m_data{function_pointer},
        m_erased_function_pointer{[](void* function_pointer, Args... args) {
          return invoke(static_cast<R (*)(Args...)>(function_pointer),
                        std::forward<Args>(args)...);
        }} {}

  template <typename F,
            typename = std::enable_if_t<!std::is_same<F, function_ref>::value>,
            typename = std::enable_if_t<is_invocable<F, Args...>::value>>
  function_ref(F&& function_object)
      : m_data{std::addressof(function_object)},
        m_erased_function_pointer{[](void* function, Args... args) -> R {
          return invoke(*static_cast<F*>(function),
                        std::forward<Args>(args)...);
        }} {}

  R operator()(Args... args) const {
    return invoke(m_erased_function_pointer, m_data,
                  std::forward<Args>(args)...);
  }

private:
  using function_pointer = R (*)(Args...);
  using erased_pointer = R (*)(void*, Args...);

  void* m_data{nullptr};
  erased_pointer m_erased_function_pointer;
};

} // namespace v1
} // namespace fub

#endif