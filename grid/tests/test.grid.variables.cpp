#include "fub/variables.hpp"

using namespace fub;

// clang-format off
struct Density {};
constexpr bool operator==(const Density&, const Density&) noexcept { return true; }
constexpr bool operator!=(const Density&, const Density&) noexcept { return false; }

static_assert(std::is_nothrow_copy_constructible<Density>::value, "");
static_assert(is_nothrow_equality_comparable<Density>::value, "");

struct Momentum {};
constexpr bool operator==(const Momentum&, const Momentum&) noexcept { return true; }
constexpr bool operator!=(const Momentum&, const Momentum&) noexcept { return false; }

static_assert(std::is_nothrow_copy_constructible<Momentum>::value, "");
static_assert(is_nothrow_equality_comparable<Momentum>::value, "");

struct Energy {};
constexpr bool operator==(const Energy&, const Energy&) noexcept { return true; }
constexpr bool operator!=(const Energy&, const Energy&) noexcept { return false; }

static_assert(std::is_nothrow_copy_constructible<Energy>::value, "");
static_assert(is_nothrow_equality_comparable<Energy>::value, "");

struct Pressure {};
constexpr bool operator==(const Pressure&, const Pressure&) noexcept { return true; }
constexpr bool operator!=(const Pressure&, const Pressure&) noexcept { return false; }

static_assert(std::is_nothrow_copy_constructible<Pressure>::value, "");
static_assert(is_nothrow_equality_comparable<Pressure>::value, "");
// clang-format on

struct State : static_collection<Density, Momentum, Energy> {};

static_assert(std::is_nothrow_copy_constructible<State>::value, "");
static_assert(is_equality_comparable<State>::value, "");
static_assert(is_nothrow_equality_comparable<State>::value, "");

static_assert(collection_traits<State>::variables(State()) ==
                  hana::make_tuple(Density{}, Momentum{}, Energy{}),
              "");

static_assert(collection_traits<Density>::variables(Density()) ==
                  hana::make_tuple(Density{}),
              "");

static_assert(collection_traits<Pressure>::variables(Pressure()) ==
                  hana::make_tuple(Pressure{}),
              "");

static_assert(collection_traits<Momentum>::variables(Momentum()) ==
                  hana::make_tuple(Momentum{}),
              "");

int main() {
  constexpr dynamic_collection collection(10);
  constexpr static_collection<Density, dynamic_collection> sc(Density{},
                                                              collection);
  static_assert(collection_traits<decltype(sc)>::size(sc) == 11, "");
}