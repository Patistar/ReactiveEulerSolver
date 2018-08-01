#ifndef FUB_CORE_NATIVE_ADDRESS_HPP
#define FUB_CORE_NATIVE_ADDRESS_HPP

namespace fub {
inline namespace v1 {

struct native_address {
  using difference_type = std::ptrdiff_t;

  void* address{nullptr};

  constexpr native_address() = default;

  template <typename T>
  constexpr explicit native_address(T* ptr)
      : address{static_cast<void*>(ptr)} {}

  template <typename T>
  constexpr std::ptrdiff_t distance_as(native_address a) const noexcept {
    return static_cast<T*>(a.address) - static_cast<T*>(address);
  }

  template <typename T>
  constexpr native_address& advance_as(std::ptrdiff_t n) noexcept {
    address = static_cast<void*>(static_cast<T*>(address) + n);
    return *this;
  }

  constexpr explicit operator void*() const noexcept { return address; }
  constexpr explicit operator const void*() const noexcept { return address; }
};

} // namespace v1
} // namespace fub

#endif
