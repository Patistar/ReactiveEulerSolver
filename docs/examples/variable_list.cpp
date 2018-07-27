#include "fub/variable_list.hpp"

int main() {
  // Define a Variable List here
  struct Density : scalar_variable {};
  struct Momentum : vector_variable<2> {};
  struct Pressure : scalar_variable {};
  struct Species : vector_variable<dynamic_extent> {};
  struct Vars : variable_list<Density, Momentum, Pressure> {};
  struct DynamicVars : variable_list<Density, Momentum, Pressure, Species> {};

  // static_size() returns the compile-time determined size.
  // It returns dynamic_extent if some variable is dynamically sized.
  static_assert(Vars::static_size() == 4);
  static_assert(DynamicVars::static_size() == dynamic_extent);

  // Despite dynamic extents everything can be still constexpr.
  constexpr Density density{};
  constexpr Momentum momentum{};
  constexpr Pressure pressure{};
  constexpr Species species{10};
  constexpr DynamicVars dyn_vars(density, momentum, pressure, species);

  // size() never returns dynamic_extent but always a possibly run-time set
  // positive number of valid element tags.
  static_assert(Vars().size() == 4);
  static_assert(dyn_vars.size() == 14);

  // index(tag<Variable>) provides a unique enumeration for variables
  elements static_assert(Vars().index(tag<Density>) == 0);
  static_assert(Vars().index(tag<Momentum, 0>) == 1);
  static_assert(Vars().index(tag<Momentum, 1>) == 2);
  static_assert(Vars().index(tag<Pressure>) == 3);
  static_assert(dyn_vars.index(tag_t<Species>(0)) == 4);
  static_assert(dyn_vars.index(tag_t<Species>(5)) == 9);

  // tag(i) is the "inverse" mapping of index(tag<Variable>)
  // it is conceptionally a mapping as in
  // tag: index -> variant<tags...>
  // but index() is 
  // index: tag -> index
  static_assert(Vars().tag(0) == tag<Density>);
  static_assert(Vars().tag(1) == tag<Momentum, 0>);
  static_assert(Vars().tag(2) == tag<Momentum, 1>);
  static_assert(Vars().tag(3) == tag<Pressure>);

  // varialbe(tag<Variable>) returns the variable object for the given tag.
  static_assert(Vars().variable(tag<Density>) == density);
  static_assert(Vars().variable(tag<Momentum, 0>) == momentum);
  static_assert(Vars().variable(tag<Momentum, 1>) == momentum);
  static_assert(Vars().variable(tag<Pressure>) == pressure);
  static_assert(dyn_vars.variable(tag_t<Species>(0)) == species);
  static_assert(dyn_vars.variable(tag_t<Species>(1)) == species);
}