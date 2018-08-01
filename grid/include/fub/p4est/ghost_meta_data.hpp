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

#ifndef FUB_P4EST_GHOST_META_DATA_HPP
#define FUB_P4EST_GHOST_META_DATA_HPP

extern "C" {
#include <p4est_ghost.h>
}

#include <memory>

namespace fub {
inline namespace v1 {
namespace p4est {

template <int Rank> struct ghost_meta_data;

template <> struct ghost_meta_data<2> {
public:
  ghost_meta_data() = default;

  /// Takes the ownership of ghost.
  explicit ghost_meta_data(p4est_t* forest)
      : m_handle{p4est_ghost_new(forest, P4EST_CONNECT_FACE)} {}

  p4est_ghost_t* operator->() const noexcept {
    return m_handle.get();
  }

  operator p4est_ghost_t*() const noexcept {
    return m_handle.get();
  }
  // operator const p4est_ghost_t*() const noexcept;

private:
  struct destroyer {
    void operator()(p4est_ghost_t* p) const noexcept {
      if (p) {
        p4est_ghost_destroy(p);
      }
    }
  };
  std::unique_ptr<p4est_ghost_t, destroyer> m_handle{nullptr};
};

} // namespace p4est
} // namespace v1
} // namespace fub

#endif