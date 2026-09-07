# The C ABI shared library (declared in morton_c.h, defined in morton_c.cpp). The Python package
# loads it via ctypes; other languages can link it directly.
add_library(peclet_morton_c SHARED morton_c.cpp)
target_link_libraries(peclet_morton_c PRIVATE morton)
target_include_directories(peclet_morton_c PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(peclet_morton_c PROPERTIES
    CXX_STANDARD 17
    POSITION_INDEPENDENT_CODE ON
    CXX_VISIBILITY_PRESET hidden)
# PECLET_MORTON_C_BUILD switches the header's export macro to dllexport on Windows.
target_compile_definitions(peclet_morton_c PRIVATE PECLET_MORTON_C_BUILD=1)
target_compile_options(peclet_morton_c PRIVATE $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-O3>)

# Drop the shared library next to the Python package so the ctypes wrapper
# finds it without an install step (LIBRARY = .so/.dylib, RUNTIME = .dll).
set_target_properties(peclet_morton_c PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python/peclet/morton
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python/peclet/morton)
