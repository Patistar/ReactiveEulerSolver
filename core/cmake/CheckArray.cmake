set(CMAKE_REQUIRED_FLAGS "-std=c++${FUB_CORE_CXX_STANDARD}")
check_cxx_source_compiles("
#include <array>

int main() {
  constexpr std::array<int, 3> array{1, 2, 3};
  constexpr int i1 = array[1];
}
"
        FUB_CORE_USE_STD_ARRAY)
