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

#ifndef FUB_READY_FUTURE_HPP
#define FUB_READY_FUTURE_HPP

#include "fub/optional.hpp"
#include "fub/type_traits.hpp"

namespace fub {
template <typename T> class ready_future {
public:
  ready_future() noexcept : m_value{} {}
  ready_future(const ready_future&) = delete;
  ready_future& operator=(const ready_future&) = delete;
  ready_future(ready_future&&) = default;
  ready_future& operator=(ready_future&&) = default;

  explicit ready_future(T&& value) : m_value{std::move(value)} {}

  explicit ready_future(const T& value) : m_value{value} {}

  ready_future(ready_future<ready_future> result)
      : ready_future([&] {
          if (result.m_value) {
            return ready_future(std::move(*result.m_value));
          }
          return ready_future();
        }()) {}

  ready_future& operator=(ready_future<ready_future> result) {
    if (result.m_value) {
      m_value.reset(std::move(*result.m_value));
    } else {
      m_value.reset();
    }
    return *this;
  }

  T get() {
    assert(m_value);
    T result{*std::move(m_value)};
    return result;
  }

  void wait() const noexcept {}

  template <typename F>
  std::enable_if_t<is_invocable<F, ready_future&&>::value,
                   ready_future<invoke_result_t<F, ready_future&&>>>
  then(F continuation) {
    assert(m_value);
    return ready_future<invoke_result_t<F, ready_future&&>>(
        fub::invoke(continuation, std::move(*this)));
  }

private:
  template <typename S> friend class ready_future;
  optional<T> m_value;
};

template <> class ready_future<void> {
public:
  ready_future() = default;
  ready_future(ready_future<ready_future>) {}
  ready_future& operator=(ready_future<ready_future>) { return *this; }

  void get() {}
  void wait() const noexcept {}

  template <typename F>
  std::enable_if_t<is_invocable<F, ready_future&&>::value,
                   ready_future<invoke_result_t<F, ready_future&&>>>
  then(F continuation) {
    return ready_future<invoke_result_t<F, ready_future&&>>(
        fub::invoke(continuation, std::move(*this)));
  }
};
} // namespace fub

#endif