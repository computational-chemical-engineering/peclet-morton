# CUDA backend for morton-arithmetic.
#
# Standalone build for a normal CUDA Toolkit install (nvcc):
#   cmake -S cuda -B cuda/build -DCMAKE_CUDA_ARCHITECTURES=90
#   cmake --build cuda/build -j
#   ctest --test-dir cuda/build --output-on-failure
#
# (Use the architecture(s) of your GPU; "90" PTX JITs onto newer cards such as
#  sm_120 / RTX 50-series. For a clang-based build without a full nvcc, see
#  build_clang.sh.)

cmake_minimum_required(VERSION 3.18)
project(morton_cuda LANGUAGES CXX CUDA)

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 90)
endif()
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# The morton core headers (sibling library) and this backend's headers.
set(MORTON_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../include" CACHE PATH
    "Path to the morton core include/ directory")

add_library(morton_cuda INTERFACE)
target_include_directories(morton_cuda INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include ${MORTON_INCLUDE_DIR})
target_compile_features(morton_cuda INTERFACE cxx_std_17)

option(MORTON_CUDA_BUILD_TESTS "Build CUDA tests/benchmarks" ON)
if(MORTON_CUDA_BUILD_TESTS AND CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    enable_testing()

    add_executable(morton_cuda_tests tests/test_cuda.cu)
    target_link_libraries(morton_cuda_tests PRIVATE morton_cuda)
    add_test(NAME morton_cuda_tests COMMAND morton_cuda_tests)

    add_executable(morton_cuda_bench bench/bench_cuda.cu)
    target_link_libraries(morton_cuda_bench PRIVATE morton_cuda)
    # Host code uses BMI2 for the CPU reference; the device pass ignores it
    # (the morton header guards PDEP behind !defined(__CUDA_ARCH__)).
    target_compile_options(morton_cuda_bench PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:-O3>
        $<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<CUDA_COMPILER_ID:NVIDIA>>:-Xcompiler=-mbmi2>)
endif()
