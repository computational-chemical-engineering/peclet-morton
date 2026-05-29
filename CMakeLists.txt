cmake_minimum_required(VERSION 3.16)
project(morton_arithmetic VERSION 0.1.0 LANGUAGES CXX)

# ---------------------------------------------------------------------------
# Header-only interface library
# ---------------------------------------------------------------------------
add_library(morton INTERFACE)
add_library(morton::morton ALIAS morton)
target_include_directories(morton INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_compile_features(morton INTERFACE cxx_std_17)

option(MORTON_ENABLE_BMI2 "Use BMI2 PDEP/PEXT for encode/decode" ON)
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mbmi2" MORTON_HAS_MBMI2)
if(MORTON_ENABLE_BMI2 AND MORTON_HAS_MBMI2)
    target_compile_options(morton INTERFACE -mbmi2)
endif()

# Only build the extras when this is the top-level project.
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    option(MORTON_BUILD_TESTS "Build the test suite" ON)
    option(MORTON_BUILD_BENCHMARKS "Build the benchmarks" ON)
    option(MORTON_BUILD_BINDINGS "Build the C ABI shared library for Python" ON)

    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
    endif()

    if(MORTON_BUILD_TESTS)
        enable_testing()
        add_subdirectory(tests)
    endif()
    if(MORTON_BUILD_BENCHMARKS)
        add_subdirectory(benchmarks)
    endif()
    if(MORTON_BUILD_BINDINGS)
        add_subdirectory(bindings)
    endif()
endif()

# ---------------------------------------------------------------------------
# Install (header-only)
# ---------------------------------------------------------------------------
include(GNUInstallDirs)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS morton EXPORT mortonTargets)
install(EXPORT mortonTargets
    FILE mortonTargets.cmake
    NAMESPACE morton::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/morton)
