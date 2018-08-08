#include "fub/variable_list.hpp"

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

// We define some variables here
// These can be used to reference to some mapped data

struct Density : fub::scalar_variable {};

template <int Rank> struct Momentum : fub::vector_variable<Rank> {};

struct Energy : fub::scalar_variable {};

struct Pressure : fub::scalar_variable {};

struct Species : fub::vector_variable<fub::dynamic_extent> {
  using fub::vector_variable<fub::dynamic_extent>::vector_variable;
};

// Here we bundle them to one variable list
// These can share a common data pointer which can be accessed into by an offset

static constexpr int Rank = 3;

struct StateVariables
    : fub::variable_list<Density, Momentum<Rank>, Energy, Species> {
  using fub::variable_list<Density, Momentum<Rank>, Energy,
                           Species>::variable_list;
};

static_assert(StateVariables::static_size() == fub::dynamic_extent, "");
static_assert(Density::static_size() == 1, "");
static_assert(Momentum<Rank>::static_size() == Rank, "");
static_assert(Energy::static_size() == 1, "");
static_assert(Species::static_size() == fub::dynamic_extent, "");

TEST_CASE("mixed variable list use cases") {
  std::ptrdiff_t N = 50;

  Density density{};
  Momentum<Rank> momentum{};
  Energy energy{};
  Species species{N};

  // Test each single variable size
  REQUIRE(density.size() == 1);
  REQUIRE(momentum.size() == Rank);
  REQUIRE(energy.size() == 1);
  REQUIRE(species.size() == N); // This is a run-time condition

  StateVariables state_variables{density, momentum, energy, species};

  SECTION("Manually test all entries") {
    // Test for correct total runtime size
    std::ptrdiff_t sum_of_sizes =
        density.size() + momentum.size() + energy.size() + species.size();
    REQUIRE(state_variables.size() == sum_of_sizes);

    // Manually test indices for each variable tag
    REQUIRE(state_variables.index(fub::tag<Density>) == 0);
    REQUIRE(state_variables.index(fub::tag<Momentum<Rank>, 0>) == 1);
    REQUIRE(state_variables.index(fub::tag<Momentum<Rank>, 1>) == 2);
    REQUIRE(state_variables.index(fub::tag<Momentum<Rank>, 2>) == 3);
    REQUIRE(state_variables.index(fub::tag<Energy>) == 4);
    // Can not find pressure
    REQUIRE(!bool(state_variables.index(fub::tag<Pressure>)));

    // Iterate through all species
    for (int i = 0; i < N; ++i) {
      REQUIRE(state_variables.index(fub::tag_t<Species>(i)) == 5 + i);
    }

    // Can not find invalid species
    REQUIRE(!bool(state_variables.index(fub::tag_t<Species>(-1))));
    REQUIRE(!bool(state_variables.index(fub::tag_t<Species>(51))));

    // variable() is the inverse of index()

    REQUIRE(state_variables.tag(0) == fub::tag<Density>);
    REQUIRE(state_variables.tag(1) == fub::tag<Momentum<Rank>, 0>);
    REQUIRE(state_variables.tag(2) == fub::tag<Momentum<Rank>, 1>);
    REQUIRE(state_variables.tag(3) == fub::tag<Momentum<Rank>, 2>);
    REQUIRE(state_variables.tag(4) == fub::tag<Energy>);

    // Iterate through all species
    for (int i = 0; i < N; ++i) {
      REQUIRE(state_variables.tag(5 + i) == fub::tag_t<Species>(i));
    }
  }

  SECTION("traverse through all state variables") {
    std::ptrdiff_t counter = 0;
    for_each(state_variables, [&](auto tag) {
      auto index = state_variables.index(tag);
      REQUIRE(index);
      auto mapped_tag = state_variables.tag(*index);
      REQUIRE(mapped_tag);
      REQUIRE(mapped_tag == tag);
      counter += 1;
    });
    REQUIRE(counter == state_variables.size());
  }
}

struct Mask : fub::scalar_variable {};
struct MaskVariables : fub::variable_list<Mask> {
  using fub::variable_list<Mask>::variable_list;
};

TEST_CASE("variable map use case") {
  constexpr int Rank = 3;
  constexpr std::ptrdiff_t N = 50;

  constexpr Density density{};
  constexpr Momentum<Rank> momentum{};
  constexpr Energy energy{};
  constexpr Species species(N);
  constexpr StateVariables state_variables(density, momentum, energy, species);

  constexpr Mask mask{};
  constexpr MaskVariables mask_variables{mask};
}