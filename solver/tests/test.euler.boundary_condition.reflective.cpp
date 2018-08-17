#include "fub/euler/boundary_condition/reflective.hpp"
#include "fub/euler/perfect_gas.hpp"
#include "fub/p4est/grid.hpp"

#include <mpi.h>

void test_reflective_boundary() {
  using namespace fub;
  constexpr int Rank = 2;
  auto grid = make_simple_grid(euler::perfect_gas_complete(),
                               dynamic_extents_t<Rank>(8, 8));
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  final_act _{[] { MPI_Finalize(); }};
  test_reflective_boundary();
}