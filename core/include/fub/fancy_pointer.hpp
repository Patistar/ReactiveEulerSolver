#ifndef FUB_CORE_SYNTHETIC_POINTER_HPP
#define FUB_CORE_SYNTHETIC_POINTER_HPP

#include "fub/type_traits.hpp"
#include "fub/utility.hpp"

#include <cassert>
#include <iterator>

namespace fub {
inline namespace v1 {

template <typename T, typename AddressModel> class fancy_pointer {
public:
  // Define these types to fulfill RandomAccessIterator concept

  using difference_type = typename AddressModel::difference_type;
  using value_type = std::remove_cv_t<T>;
  using reference = std::add_lvalue_reference_t<T>;
  using pointer = T*;
  using iterator_category = std::random_access_iterator_tag;

  // Constructors and Assignments to fulfill the NullablePointer concept

  constexpr fancy_pointer() = default;

  constexpr fancy_pointer(const fancy_pointer& p) noexcept = default;

  constexpr fancy_pointer(fancy_pointer&& p) noexcept = default;

  constexpr fancy_pointer(std::nullptr_t) noexcept // NOLINT
      : fancy_pointer() {}

  constexpr fancy_pointer& operator=(std::nullptr_t) noexcept {
    m_address = AddressModel{};
    return *this;
  }

  constexpr fancy_pointer(T* native) : m_address(native) {}

  template <typename... Args, typename std::enable_if_t<std::is_constructible<
                                  AddressModel, Args...>::value>>
  constexpr fancy_pointer(Args&&... args)
      : m_address(std::forward<Args>(args)...) {}

  template <typename S,
            typename = std::enable_if_t<std::is_convertible<S*, T*>::value>>
  constexpr fancy_pointer(const fancy_pointer<S, AddressModel>& p)
      : m_address{p.m_address} {}

  // Contextually Convertible to bool

  explicit constexpr operator bool() const noexcept {
    return static_cast<const void*>(m_address) != nullptr;
  }

  // Access Operators

  constexpr pointer operator->() const noexcept {
    assert(static_cast<const void*>(m_address) != nullptr);
    return static_cast<T*>(static_cast<void*>(m_address));
  }

  constexpr reference operator*() const noexcept {
    return *static_cast<T*>(static_cast<void*>(m_address));
  }

  constexpr reference operator[](difference_type n) const noexcept {
    AddressModel address = m_address;
    address.template advance_as<T>(n);
    return *static_cast<T*>(static_cast<void*>(address));
  }

  // Operators for the RandomAccessIterator

  constexpr fancy_pointer& operator+=(difference_type n) noexcept {
    m_address.template advance_as<T>(n);
    return *this;
  }
  constexpr fancy_pointer& operator-=(difference_type n) noexcept {
    m_address.template advance_as<T>(-n);
    return *this;
  }

  friend constexpr difference_type operator-(fancy_pointer a,
                                             fancy_pointer b) noexcept {
    return b.m_address.template distance_as<T>(a.m_address);
  }

  constexpr explicit operator T* () const noexcept {
    return static_cast<T*>(static_cast<void*>(m_address));
  }

  // constexpr operator void* () const noexcept {
  //   return static_cast<void*>(m_address);
  // }

private:
  AddressModel m_address{};
};

template <typename AddressModel> class fancy_pointer<void, AddressModel> {
public:
  // Define these types to fulfill RandomAccessIterator concept

  using difference_type = typename AddressModel::difference_type;
  using value_type = void;
  using reference = void;
  using pointer = void*;
  using iterator_category = std::random_access_iterator_tag;

  // Constructors and Assignments to fulfill the NullablePointer concept

  constexpr fancy_pointer() = default;

  constexpr fancy_pointer(std::nullptr_t) noexcept // NOLINT
      : fancy_pointer() {}

  constexpr fancy_pointer& operator=(std::nullptr_t) noexcept {
    m_address = AddressModel{};
    return *this;
  }

  template <typename... Args, typename std::enable_if_t<std::is_constructible<
                                  AddressModel, Args...>::value>>
  constexpr explicit fancy_pointer(Args&&... args) noexcept
      : m_address{std::forward<Args>(args)...} {}

  constexpr fancy_pointer(const fancy_pointer& p) noexcept = default;

  template <typename S>
  constexpr explicit fancy_pointer(
      const fancy_pointer<S, AddressModel>& p) noexcept
      : m_address{p.m_address} {}

  // Contextually Convertible to bool

  explicit constexpr operator bool() const noexcept {
    return m_address.is_null();
  }

  // Access Operators

  constexpr pointer operator->() const noexcept {
    assert(!m_address.is_null());
    return static_cast<pointer>(static_cast<void*>(m_address));
  }

  // Operators for the RandomAccessIterator

  constexpr fancy_pointer& operator+=(difference_type n) noexcept {
    m_address.advance(n);
    return *this;
  }
  constexpr fancy_pointer& operator-=(difference_type n) noexcept {
    m_address.advance(-n);
    return *this;
  }

  friend constexpr difference_type operator-(fancy_pointer a,
                                             fancy_pointer b) noexcept {
    return a.m_address.distance(b.m_address);
  }

  constexpr operator void*() const noexcept {
    return static_cast<void*>(m_address);
  }

private:
  AddressModel m_address{};
};

template <typename T, typename AM>
constexpr fancy_pointer<T, AM>
operator+(fancy_pointer<T, AM> p,
          typename fancy_pointer<T, AM>::difference_type n) noexcept {
  fancy_pointer<T, AM> temp = p;
  return temp += n;
}

template <typename T, typename AM>
constexpr fancy_pointer<T, AM>
operator+(typename fancy_pointer<T, AM>::difference_type n,
          fancy_pointer<T, AM> p) noexcept {
  return p + n;
}

template <typename T, typename AM>
constexpr fancy_pointer<T, AM>
operator-(fancy_pointer<T, AM> p,
          typename fancy_pointer<T, AM>::difference_type n) noexcept {
  fancy_pointer<T, AM> temp = p;
  return temp -= n;
}

template <typename T, typename AM>
constexpr bool operator==(fancy_pointer<T, AM> p, std::nullptr_t) noexcept {
  return static_cast<bool>(p);
}

template <typename T, typename AM>
constexpr bool operator==(std::nullptr_t, fancy_pointer<T, AM> p) noexcept {
  return p == nullptr;
}

template <typename T, typename AM>
constexpr bool operator!=(fancy_pointer<T, AM> p, std::nullptr_t) noexcept {
  return !(p == nullptr);
}

template <typename T, typename AM>
constexpr bool operator!=(std::nullptr_t, fancy_pointer<T, AM> p) noexcept {
  return !(nullptr == p);
}

template <typename T, typename AM>
constexpr bool operator==(fancy_pointer<T, AM> a,
                          nodeduce_t<fancy_pointer<T, AM>> b) noexcept {
  return typename fancy_pointer<T, AM>::difference_type{0} == (b - a);
}

template <typename T, typename AM>
constexpr bool operator!=(fancy_pointer<T, AM> a,
                          nodeduce_t<fancy_pointer<T, AM>> b) noexcept {
  return !(a == b);
}

template <typename T, typename AM>
constexpr bool operator<(fancy_pointer<T, AM> a,
                         nodeduce_t<fancy_pointer<T, AM>> b) noexcept {
  return typename fancy_pointer<T, AM>::difference_type{0} < (b - a);
}

template <typename T, typename AM>
constexpr bool operator<=(fancy_pointer<T, AM> a,
                          nodeduce_t<fancy_pointer<T, AM>> b) noexcept {
  return typename fancy_pointer<T, AM>::difference_type{0} <= (b - a);
}

template <typename T, typename AM>
constexpr bool operator>(fancy_pointer<T, AM> a,
                         nodeduce_t<fancy_pointer<T, AM>> b) noexcept {
  return typename fancy_pointer<T, AM>::difference_type{0} > (b - a);
}

template <typename T, typename AM>
constexpr bool operator>=(fancy_pointer<T, AM> a,
                          nodeduce_t<fancy_pointer<T, AM>> b) noexcept {
  return typename fancy_pointer<T, AM>::difference_type{0} >= (b - a);
}

template <typename T, typename AM>
constexpr fancy_pointer<T, AM>& operator++(fancy_pointer<T, AM>& p) noexcept {
  return p += 1;
}

template <typename T, typename AM>
constexpr fancy_pointer<T, AM>& operator--(fancy_pointer<T, AM>& p) noexcept {
  return p -= 1;
}

template <typename T, typename AM>
constexpr fancy_pointer<T, AM> operator++(fancy_pointer<T, AM> p,
                                          int) noexcept {
  fancy_pointer<T, AM> q = p;
  return q += 1;
}

template <typename T, typename AM>
constexpr fancy_pointer<T, AM> operator--(fancy_pointer<T, AM> p,
                                          int) noexcept {
  fancy_pointer<T, AM> q = p;
  return q -= 1;
}

} // namespace v1
} // namespace fub

#endif
