// morton_cuda/morton_cuda.cuh
//
// CUDA backend for the morton-arithmetic library. The per-element work is the
// ordinary Morton<Dim,Bits> code path (made __host__ __device__ by MORTON_HD in
// morton/morton.hpp), so the GPU kernels share the exact, already-tested logic
// of the CPU library -- there is no second implementation to keep in sync.
//
// Two layers:
//   *_device(...)  operate on device pointers the caller owns (compose into
//                  larger pipelines, pick streams).
//   *_host(...)    convenience wrappers that allocate, copy H<->D, launch and
//                  free -- handy for one-shot calls and tests.
//
// Encode/decode on the GPU necessarily use the software (PDEP-free) bit path;
// there is no PDEP/PEXT on NVIDIA hardware.
//
// SPDX-License-Identifier: MIT

#ifndef MORTON_CUDA_MORTON_CUDA_CUH
#define MORTON_CUDA_MORTON_CUDA_CUH

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>

#include "morton/morton.hpp"

namespace morton {
namespace cuda {

inline void check(cudaError_t e, const char* what) {
    if (e != cudaSuccess)
        throw std::runtime_error(std::string("morton::cuda: ") + what + ": " +
                                 cudaGetErrorString(e));
}

namespace detail {
constexpr int block = 256;
inline int grid(std::size_t n) { return int((n + block - 1) / block); }
}  // namespace detail

// ---------------------------------------------------------------------------
// Kernels
// ---------------------------------------------------------------------------
template <unsigned Bits>
__global__ void encode2_kernel(const typename Morton<2, Bits>::coord_type* x,
                               const typename Morton<2, Bits>::coord_type* y,
                               typename Morton<2, Bits>::code_type* out, std::size_t n) {
    std::size_t i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < n) out[i] = Morton<2, Bits>::encode(x[i], y[i]).code();
}

template <unsigned Bits>
__global__ void decode2_kernel(const typename Morton<2, Bits>::code_type* code,
                               typename Morton<2, Bits>::coord_type* x,
                               typename Morton<2, Bits>::coord_type* y, std::size_t n) {
    std::size_t i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < n) {
        Morton<2, Bits> m = Morton<2, Bits>::from_code(code[i]);
        x[i] = m.get(0);
        y[i] = m.get(1);
    }
}

template <unsigned Bits>
__global__ void encode3_kernel(const typename Morton<3, Bits>::coord_type* x,
                               const typename Morton<3, Bits>::coord_type* y,
                               const typename Morton<3, Bits>::coord_type* z,
                               typename Morton<3, Bits>::code_type* out, std::size_t n) {
    std::size_t i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < n) out[i] = Morton<3, Bits>::encode(x[i], y[i], z[i]).code();
}

template <unsigned Bits>
__global__ void decode3_kernel(const typename Morton<3, Bits>::code_type* code,
                               typename Morton<3, Bits>::coord_type* x,
                               typename Morton<3, Bits>::coord_type* y,
                               typename Morton<3, Bits>::coord_type* z, std::size_t n) {
    std::size_t i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < n) {
        Morton<3, Bits> m = Morton<3, Bits>::from_code(code[i]);
        x[i] = m.get(0);
        y[i] = m.get(1);
        z[i] = m.get(2);
    }
}

// Add a signed delta to one axis of each code (wraps mod 2^Bits).
template <unsigned Dim, unsigned Bits>
__global__ void add_kernel(const typename Morton<Dim, Bits>::code_type* in,
                           typename Morton<Dim, Bits>::code_type* out, std::size_t n,
                           unsigned axis, long long delta) {
    using M = Morton<Dim, Bits>;
    using C = typename M::coord_type;
    std::size_t i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < n) {
        M m = M::from_code(in[i]);
        if (delta >= 0)
            m.add(axis, C((unsigned long long)delta));
        else
            m.sub(axis, C((unsigned long long)(-delta)));
        out[i] = m.code();
    }
}

// ---------------------------------------------------------------------------
// Device-pointer API (caller owns the device memory)
// ---------------------------------------------------------------------------
template <unsigned Bits>
void encode2_device(const typename Morton<2, Bits>::coord_type* dx,
                    const typename Morton<2, Bits>::coord_type* dy,
                    typename Morton<2, Bits>::code_type* dout, std::size_t n,
                    cudaStream_t s = 0) {
    if (n) encode2_kernel<Bits><<<detail::grid(n), detail::block, 0, s>>>(dx, dy, dout, n);
}
template <unsigned Bits>
void decode2_device(const typename Morton<2, Bits>::code_type* dcode,
                    typename Morton<2, Bits>::coord_type* dx,
                    typename Morton<2, Bits>::coord_type* dy, std::size_t n, cudaStream_t s = 0) {
    if (n) decode2_kernel<Bits><<<detail::grid(n), detail::block, 0, s>>>(dcode, dx, dy, n);
}
template <unsigned Bits>
void encode3_device(const typename Morton<3, Bits>::coord_type* dx,
                    const typename Morton<3, Bits>::coord_type* dy,
                    const typename Morton<3, Bits>::coord_type* dz,
                    typename Morton<3, Bits>::code_type* dout, std::size_t n, cudaStream_t s = 0) {
    if (n) encode3_kernel<Bits><<<detail::grid(n), detail::block, 0, s>>>(dx, dy, dz, dout, n);
}
template <unsigned Dim, unsigned Bits>
void add_device(const typename Morton<Dim, Bits>::code_type* din,
                typename Morton<Dim, Bits>::code_type* dout, std::size_t n, unsigned axis,
                long long delta, cudaStream_t s = 0) {
    if (n) add_kernel<Dim, Bits><<<detail::grid(n), detail::block, 0, s>>>(din, dout, n, axis, delta);
}

// ---------------------------------------------------------------------------
// Host-array convenience API (alloc + copy + launch + copy + free)
// ---------------------------------------------------------------------------
template <unsigned Bits>
void encode2_host(const typename Morton<2, Bits>::coord_type* x,
                  const typename Morton<2, Bits>::coord_type* y,
                  typename Morton<2, Bits>::code_type* out, std::size_t n) {
    using C = typename Morton<2, Bits>::coord_type;
    using K = typename Morton<2, Bits>::code_type;
    C *dx = nullptr, *dy = nullptr;
    K* dout = nullptr;
    check(cudaMalloc(&dx, n * sizeof(C)), "malloc x");
    check(cudaMalloc(&dy, n * sizeof(C)), "malloc y");
    check(cudaMalloc(&dout, n * sizeof(K)), "malloc out");
    check(cudaMemcpy(dx, x, n * sizeof(C), cudaMemcpyHostToDevice), "copy x");
    check(cudaMemcpy(dy, y, n * sizeof(C), cudaMemcpyHostToDevice), "copy y");
    encode2_device<Bits>(dx, dy, dout, n);
    check(cudaGetLastError(), "encode2 launch");
    check(cudaMemcpy(out, dout, n * sizeof(K), cudaMemcpyDeviceToHost), "copy out");
    cudaFree(dx);
    cudaFree(dy);
    cudaFree(dout);
}

template <unsigned Bits>
void decode2_host(const typename Morton<2, Bits>::code_type* code,
                  typename Morton<2, Bits>::coord_type* x,
                  typename Morton<2, Bits>::coord_type* y, std::size_t n) {
    using C = typename Morton<2, Bits>::coord_type;
    using K = typename Morton<2, Bits>::code_type;
    K* dcode = nullptr;
    C *dx = nullptr, *dy = nullptr;
    check(cudaMalloc(&dcode, n * sizeof(K)), "malloc code");
    check(cudaMalloc(&dx, n * sizeof(C)), "malloc x");
    check(cudaMalloc(&dy, n * sizeof(C)), "malloc y");
    check(cudaMemcpy(dcode, code, n * sizeof(K), cudaMemcpyHostToDevice), "copy code");
    decode2_device<Bits>(dcode, dx, dy, n);
    check(cudaGetLastError(), "decode2 launch");
    check(cudaMemcpy(x, dx, n * sizeof(C), cudaMemcpyDeviceToHost), "copy x");
    check(cudaMemcpy(y, dy, n * sizeof(C), cudaMemcpyDeviceToHost), "copy y");
    cudaFree(dcode);
    cudaFree(dx);
    cudaFree(dy);
}

template <unsigned Bits>
void encode3_host(const typename Morton<3, Bits>::coord_type* x,
                  const typename Morton<3, Bits>::coord_type* y,
                  const typename Morton<3, Bits>::coord_type* z,
                  typename Morton<3, Bits>::code_type* out, std::size_t n) {
    using C = typename Morton<3, Bits>::coord_type;
    using K = typename Morton<3, Bits>::code_type;
    C *dx = nullptr, *dy = nullptr, *dz = nullptr;
    K* dout = nullptr;
    check(cudaMalloc(&dx, n * sizeof(C)), "malloc x");
    check(cudaMalloc(&dy, n * sizeof(C)), "malloc y");
    check(cudaMalloc(&dz, n * sizeof(C)), "malloc z");
    check(cudaMalloc(&dout, n * sizeof(K)), "malloc out");
    check(cudaMemcpy(dx, x, n * sizeof(C), cudaMemcpyHostToDevice), "copy x");
    check(cudaMemcpy(dy, y, n * sizeof(C), cudaMemcpyHostToDevice), "copy y");
    check(cudaMemcpy(dz, z, n * sizeof(C), cudaMemcpyHostToDevice), "copy z");
    encode3_device<Bits>(dx, dy, dz, dout, n);
    check(cudaGetLastError(), "encode3 launch");
    check(cudaMemcpy(out, dout, n * sizeof(K), cudaMemcpyDeviceToHost), "copy out");
    cudaFree(dx);
    cudaFree(dy);
    cudaFree(dz);
    cudaFree(dout);
}

template <unsigned Dim, unsigned Bits>
void add_host(const typename Morton<Dim, Bits>::code_type* in,
              typename Morton<Dim, Bits>::code_type* out, std::size_t n, unsigned axis,
              long long delta) {
    using K = typename Morton<Dim, Bits>::code_type;
    K *din = nullptr, *dout = nullptr;
    check(cudaMalloc(&din, n * sizeof(K)), "malloc in");
    check(cudaMalloc(&dout, n * sizeof(K)), "malloc out");
    check(cudaMemcpy(din, in, n * sizeof(K), cudaMemcpyHostToDevice), "copy in");
    add_device<Dim, Bits>(din, dout, n, axis, delta);
    check(cudaGetLastError(), "add launch");
    check(cudaMemcpy(out, dout, n * sizeof(K), cudaMemcpyDeviceToHost), "copy out");
    cudaFree(din);
    cudaFree(dout);
}

}  // namespace cuda
}  // namespace morton

#endif  // MORTON_CUDA_MORTON_CUDA_CUH
