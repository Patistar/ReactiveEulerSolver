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

#ifndef FUB_INTERVALMAP_HPP
#define FUB_INTERVALMAP_HPP

#include "fub/type_traits.hpp"
#include <limits>
#include <map>

namespace fub {

/// @brief This class maps intervals of type [Key, Key) onto values of type T.
///
/// This map is initialized as to map all keys to a default value of T, which
/// might be specified in the constructor. The default value can not be changed
/// after construction.
///
/// @note Key needs to specialize std::numeric_limits<Key> and needs to be
/// LessThanComparable
template <typename Key, typename T, typename MapT = std::map<Key, T>>
class interval_map {
private:
  MapT m_map;

public:
  using allocator = typename MapT::allocator_type;
  using key_compare = typename MapT::key_compare;

  explicit interval_map(T default_value = T{})
      : m_map{{std::make_pair(std::numeric_limits<Key>::lowest(),
                              std::move(default_value))}} {}

  // DEFAULT CONSTRUCTORS

  interval_map(const interval_map&) = default;
  interval_map& operator=(const interval_map&) = default;

  interval_map(interval_map&&) = default;
  interval_map& operator=(interval_map&&) = default;

  ~interval_map() = default;

  // MORE CONSTRUCTORS

  explicit interval_map(const key_compare& comp,
                       const allocator& alloc = allocator()) noexcept
      : m_map{comp, alloc} {}

  interval_map(const interval_map& other, const allocator& alloc)
      : m_map{other.m_map, alloc} {}

  interval_map(interval_map&& other, const allocator& alloc) noexcept(
      std::is_nothrow_move_constructible<MapT>::value)
      : m_map{std::move(other.m_map), alloc} {}

  // OBSERVERS

  const MapT& get_map() const noexcept { return m_map; }

  const T& get_default_value() const noexcept {
    assert(m_map.size() > 0);
    return std::prev(m_map.end())->second;
  }

  // ACCESSORS

  const T& operator[](const Key& key) const noexcept {
    auto ub = m_map.upper_bound(key);
    if (ub == m_map.begin()) {
      assert(m_map.size() > 0);
      return get_default_value();
    }
    return (--ub)->second;
  }

  // CAPACITY

  bool empty() const noexcept { return m_map.empty(); }

  // MODIFIERS

  /// @brief assigns a `value` to a given interval [lower, upper).
  void insert(Key lower, Key upper, T value) {
    // If an empty interval is specified, we do nothing.
    if (!invoke(m_map.key_comp(), lower, upper)) {
      return;
    }
    auto last = m_map.upper_bound(upper);
    // this assertion holds since we have lower < upper and m_map.begin()->first
    // is the lowest key
    assert(last != m_map.begin());
    const T& value_before = std::prev(last)->second;
    last = m_map.insert(last, std::make_pair(std::move(upper), value_before));
    auto first = m_map.lower_bound(lower);
    m_map.erase(first, last);
    first =
        m_map.insert(last, std::make_pair(std::move(lower), std::move(value)));
  }
};

template <typename Key, typename T, typename Map>
bool operator==(const interval_map<Key, T, Map>& lhs,
                const interval_map<Key, T, Map>& rhs) {
  return lhs.map() == rhs.map();
}

template <typename Key, typename T, typename Map>
bool operator!=(const interval_map<Key, T, Map>& lhs,
                const interval_map<Key, T, Map>& rhs) {
  return lhs.map() != rhs.map();
}

} // namespace fub

#endif // !INTERVALMAP_HPP
