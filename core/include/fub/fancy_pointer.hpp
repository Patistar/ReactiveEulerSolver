#ifndef FUB_CORE_SYNTHETIC_POINTER_HPP
#define FUB_CORE_SYNTHETIC_POINTER_HPP

#include "fub/type_traits.hpp"
#include "fub/utility.hpp"

namespace fub {

template <typename T, typename AddressModel> class synthetic_pointer {
public:
  // Define these types to fulfill RandomAccessIterator concept

  using difference_type = typename AddressModel::difference_type;
  using value_type = std::remove_cv_t<T>;
  using reference = std::add_lvalue_reference_t<T>;
  using pointer = synthetic_pointer;
  using iterator_category = std::random_access_iterator_tag;

  // Constructors and Assignments to fulfill the NullablePointer concept

  constexpr synthetic_pointer() = default;
  constexpr synthetic_pointer(std::nullptr_t) noexcept // NOLINT
      : synthetic_pointer() {}
  constexpr synthetic_pointer& operator=(std::nullptr_t) noexcept {
    m_address = AddressModel{};
  }

  // Contextually Convertible to bool

  explicit constexpr operator bool() const noexcept {
    return m_address.is_null();
  }

  // Access Operators

  constexpr T* operator->() const noexcept {
    assert(!m_address.is_null());
    return static_cast<T*>(static_cast<void*>(m_address));
  }

  constexpr reference operator*(difference_type n) const noexcept {
    return *static_cast<T*>(static_cast<void*>(m_address));
  }

  constexpr reference operator[](difference_type n) const noexcept {
    AddressModel address = m_address;
    address.advance(n);
    return *static_cast<T*>(static_cast<void*>(address));
  }

  // Equality Operator

  friend constexpr bool operator==(synthetic_pointer p,
                                   std::nullptr_t) noexcept {
    return p.m_address.is_null();
  }
  friend constexpr bool operator==(std::nullptr_t,
                                   synthetic_pointer p) noexcept {
    return p == nullptr;
  }

  // Operators for the RandomAccessIterator

  constexpr synthetic_pointer& operator+=(difference_type n) noexcept {
    m_address.advance(n);
    return *this;
  }
  constexpr synthetic_pointer& operator-=(difference_type n) noexcept {
    m_address.advance(-n);
    return *this;
  }

  friend constexpr difference_type operator-(synthetic_pointer a,
                                             synthetic_pointer b) noexcept {
    return a.m_address.distance(b.m_address);
  }

private:
  AddressModel m_address{};
};

template <typename T, typename AM>
constexpr bool operator!=(synthetic_pointer<T, AM> p, std::nullptr_t) noexcept {
  return !(p == nullptr);
}

template <typename T, typename AM>
constexpr bool operator!=(std::nullptr_t, synthetic_pointer<T, AM> p) noexcept {
  return !(nullptr == p);
}

template <typename T, typename AM>
constexpr synthetic_pointer<T, AM>
operator+(synthetic_pointer<T, AM> p,
          typename synthetic_pointer<T, AM>::difference_type n) noexcept {
  synthetic_pointer<T, AM> temp = p;
  return temp += n;
}

template <typename T, typename AM>
constexpr synthetic_pointer<T, AM>
operator+(typename synthetic_pointer<T, AM>::difference_type n,
          synthetic_pointer<T, AM> p) noexcept {
  return p + n;
}

template <typename T, typename AM>
constexpr synthetic_pointer<T, AM>
operator-(synthetic_pointer<T, AM> p,
          typename synthetic_pointer<T, AM>::difference_type n) noexcept {
  synthetic_pointer<T, AM> temp = p;
  return temp -= n;
}

template <typename T, typename AM>
constexpr bool operator<(synthetic_pointer<T, AM> a,
                         nodeduce_t<synthetic_pointer<T, AM>> b) noexcept {
  return typename synthetic_pointer<T, AM>::difference_type{0} < (b - a);
}

template <typename T, typename AM>
constexpr bool operator<=(synthetic_pointer<T, AM> a,
                          nodeduce_t<synthetic_pointer<T, AM>> b) noexcept {
  return typename synthetic_pointer<T, AM>::difference_type{0} <= (b - a);
}

template <typename T, typename AM>
constexpr bool operator>(synthetic_pointer<T, AM> a,
                         nodeduce_t<synthetic_pointer<T, AM>> b) noexcept {
  return typename synthetic_pointer<T, AM>::difference_type{0} > (b - a);
}

template <typename T, typename AM>
constexpr bool operator>=(synthetic_pointer<T, AM> a,
                          nodeduce_t<synthetic_pointer<T, AM>> b) noexcept {
  return typename synthetic_pointer<T, AM>::difference_type{0} >= (b - a);
}

template <typename T, typename AM>
constexpr synthetic_pointer<T, AM>&
operator++(synthetic_pointer<T, AM>& p) noexcept {
  return p += 1;
}

template <typename T, typename AM>
constexpr synthetic_pointer<T, AM>&
operator--(synthetic_pointer<T, AM>& p) noexcept {
  return p -= 1;
}

template <typename T, typename AM>
constexpr synthetic_pointer<T, AM> operator++(synthetic_pointer<T, AM> p,
                                              int) noexcept {
  synthetic_pointer<T, AM> q = p;
  return q += 1;
}

template <typename T, typename AM>
constexpr synthetic_pointer<T, AM> operator--(synthetic_pointer<T, AM> p,
                                              int) noexcept {
  synthetic_pointer<T, AM> q = p;
  return q -= 1;
}

} // namespace fub

#endif
