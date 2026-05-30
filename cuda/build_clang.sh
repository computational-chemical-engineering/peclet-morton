#!/usr/bin/env bash
# Build + run the CUDA backend with clang as the CUDA compiler.
#
# Useful when you don't have a full `nvcc` (e.g. only the pip redistributable
# pieces: ptxas + libdevice + cudart) or hit nvcc/glibc header clashes. clang
# does the CUDA front-end itself and only needs ptxas/fatbinary/nvlink + headers
# in a toolkit-shaped directory pointed to by CUDA_PATH.
#
# Usage:
#   CUDA_PATH=/path/to/cuda GPU_ARCH=sm_90 ./cuda/build_clang.sh
#
# GPU_ARCH defaults to sm_90, whose PTX JITs onto newer GPUs (e.g. sm_120 /
# RTX 50-series) via the driver.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
CUDA_PATH="${CUDA_PATH:?set CUDA_PATH to a toolkit dir (bin/ptxas, nvvm/libdevice, include, lib64)}"
GPU_ARCH="${GPU_ARCH:-sm_90}"
CXX="${CXX:-clang++}"

common=(-std=c++17 -x cuda
        -I"$ROOT/include" -I"$HERE/include"
        --cuda-path="$CUDA_PATH" --cuda-gpu-arch="$GPU_ARCH" --no-cuda-version-check
        -L"$CUDA_PATH/lib64" -lcudart -ldl -lrt -pthread)

mkdir -p "$HERE/build"
echo "==> building tests"
"$CXX" -O2 "${common[@]}" "$HERE/tests/test_cuda.cu" -o "$HERE/build/morton_cuda_tests"
echo "==> building benchmark"
"$CXX" -O3 -mbmi2 "${common[@]}" "$HERE/bench/bench_cuda.cu" -o "$HERE/build/morton_cuda_bench"

echo "==> running tests"
LD_LIBRARY_PATH="$CUDA_PATH/lib64" "$HERE/build/morton_cuda_tests"
echo "==> running benchmark"
LD_LIBRARY_PATH="$CUDA_PATH/lib64" "$HERE/build/morton_cuda_bench"
