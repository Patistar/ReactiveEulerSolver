set(CMAKE_REQUIRED_FLAGS "-std=c++${FUB_CORE_CXX_STANDARD}")
check_cxx_source_compiles("
#include <functional>

int f() { return 0; }
int main() { return std::invoke(f); }"
        FUB_WITH_STD_CLAMP)
