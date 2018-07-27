set(CMAKE_REQUIRED_FLAGS "-std=c++${FUB_CORE_CXX_STANDARD}")
check_cxx_source_compiles("
#include <tuple>

int main() {
  return std::apply(
      [](auto... xs) {
        int array[]{xs...};
        return array[0];
      },
      std::make_tuple(0, 0));
}"
        FUB_WITH_STD_APPLY)

