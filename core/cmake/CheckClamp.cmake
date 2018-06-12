set(CMAKE_REQUIRED_FLAGS "-std=c++${FUB_CORE_CXX_STANDARD}")
check_cxx_source_compiles("
#include <algorithm>
int main() { return std::clamp(0, 0, 0); }"
        FUB_CORE_USE_STD_CLAMP)

