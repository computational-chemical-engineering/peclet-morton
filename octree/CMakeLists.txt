# Standalone build for the morton-octree project.
#
# It depends on the installed `morton` library:
#   cmake -S . -B build -Dmorton_DIR=<prefix>/lib/cmake/morton
# or, while both live in this monorepo, point at the sibling source tree:
#   cmake -S . -B build -DMORTON_INCLUDE_DIR=../include
#
# This is an early scaffold; see PLAN.md for where it is going.

cmake_minimum_required(VERSION 3.16)
project(morton_octree VERSION 0.0.1 LANGUAGES CXX)

add_library(morton_octree INTERFACE)
add_library(morton_octree::morton_octree ALIAS morton_octree)
target_include_directories(morton_octree INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_compile_features(morton_octree INTERFACE cxx_std_17)

# Dependency on the standalone morton library.
find_package(morton CONFIG QUIET)
if(morton_FOUND)
    target_link_libraries(morton_octree INTERFACE morton::morton)
elseif(DEFINED MORTON_INCLUDE_DIR)
    target_include_directories(morton_octree INTERFACE ${MORTON_INCLUDE_DIR})
else()
    message(FATAL_ERROR
        "morton not found: install it and pass -Dmorton_DIR=..., "
        "or pass -DMORTON_INCLUDE_DIR=<path to morton include/>.")
endif()

option(MORTON_OCTREE_BUILD_TESTS "Build octree tests" ON)
if(MORTON_OCTREE_BUILD_TESTS AND CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    # doctest single header: pass -DDOCTEST_INCLUDE_DIR=... (defaults to the
    # monorepo's vendored copy when present).
    set(DOCTEST_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../third_party"
        CACHE PATH "Directory containing doctest.h")
    enable_testing()
    add_executable(morton_octree_tests tests/test_main.cpp tests/test_octree.cpp)
    target_link_libraries(morton_octree_tests PRIVATE morton_octree)
    target_include_directories(morton_octree_tests PRIVATE ${DOCTEST_INCLUDE_DIR})
    add_test(NAME morton_octree_tests COMMAND morton_octree_tests)
endif()
