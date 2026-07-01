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
    # Only GCC/Clang understand -mbmi2; guard so the flag is safe to propagate
    # to consumers of the installed package (MSVC/other compilers ignore it).
    target_compile_options(morton INTERFACE
        $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,AppleClang>:-mbmi2>)
endif()

# Single-binary runtime ISA dispatch (x86-64 GCC/Clang). Builds *without*
# -mbmi2 but still uses PDEP/PEXT when the running CPU has BMI2 (software path
# otherwise), via per-function `target` attributes. Intended to be paired with
# MORTON_ENABLE_BMI2=OFF (it self-disables when -mbmi2 is on). The AVX-512 batch
# kernels are always runtime-dispatched and need no option. This lets a portable
# distributable binary (e.g. a wheel) keep the BMI2/AVX-512 fast paths.
option(MORTON_ENABLE_RUNTIME_DISPATCH
    "Single-binary runtime BMI2 dispatch (build without -mbmi2)" OFF)
if(MORTON_ENABLE_RUNTIME_DISPATCH)
    target_compile_definitions(morton INTERFACE MORTON_ENABLE_RUNTIME_DISPATCH=1)
endif()

# Portable Kokkos backend (include/morton/kokkos.hpp). Opt-in: the plain build
# (and the wheels) never touch Kokkos and stay bit-identical. When ON, Kokkos is
# found via find_package(CONFIG) against the suite's bootstrapped install prefix
# (extern/install/<backend>, from tools/bootstrap_deps.sh) or a cluster module --
# the same single mechanism the rest of the suite uses. The backend (CUDA / HIP /
# OpenMP / Serial) and architecture come from that Kokkos build, not from here.
option(MORTON_ENABLE_KOKKOS "Build the portable Kokkos backend (and its tests/bench)" OFF)
if(MORTON_ENABLE_KOKKOS)
    find_package(Kokkos CONFIG REQUIRED)
    message(STATUS "morton: Kokkos ${Kokkos_VERSION} (backends: ${Kokkos_DEVICES})")
    target_compile_definitions(morton INTERFACE MORTON_ENABLE_KOKKOS=1)
    target_link_libraries(morton INTERFACE Kokkos::kokkos)
endif()

# Wheel mode (driven by scikit-build-core via the root pyproject.toml): build
# only the Python extension and install it into the package directory.
option(MORTON_WHEEL "Build only the Python extension for a wheel" OFF)
if(MORTON_WHEEL)
    add_subdirectory(bindings)
    install(TARGETS mortonarith_c LIBRARY DESTINATION peclet/morton)
    return()
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
        if(MORTON_ENABLE_KOKKOS)
            add_subdirectory(tests/kokkos)
        endif()
    endif()
    if(MORTON_BUILD_BENCHMARKS)
        add_subdirectory(benchmarks)
    endif()
    if(MORTON_BUILD_BINDINGS)
        add_subdirectory(bindings)
    endif()

    # `cmake --build build --target docs` (if Doxygen is installed)
    find_package(Doxygen)
    if(DOXYGEN_FOUND)
        add_custom_target(docs
            COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_SOURCE_DIR}/docs/Doxyfile
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM)
    endif()
endif()

# ---------------------------------------------------------------------------
# Install + CMake package config so consumers can `find_package(morton)`
# ---------------------------------------------------------------------------
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(MORTON_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/morton)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS morton EXPORT mortonTargets)
install(EXPORT mortonTargets
    FILE morton-targets.cmake
    NAMESPACE morton::
    DESTINATION ${MORTON_CMAKE_DIR})

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/morton-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/morton-config.cmake
    INSTALL_DESTINATION ${MORTON_CMAKE_DIR})
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/morton-config-version.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
    ARCH_INDEPENDENT)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/morton-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/morton-config-version.cmake
    DESTINATION ${MORTON_CMAKE_DIR})
