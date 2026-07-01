add_library(mortonarith_c SHARED morton_c.cpp)
target_link_libraries(mortonarith_c PRIVATE morton)
set_target_properties(mortonarith_c PROPERTIES
    CXX_STANDARD 17
    POSITION_INDEPENDENT_CODE ON)
target_compile_options(mortonarith_c PRIVATE -O3)

# Drop the shared library next to the Python package so the ctypes wrapper
# finds it without an install step.
set_target_properties(mortonarith_c PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python/peclet/morton)
