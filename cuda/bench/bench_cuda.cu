// Throughput benchmark for the CUDA backend.
//
// Reports three numbers for 2D 32-bit encode:
//   - GPU device-only  (data already resident on the GPU; pure kernel),
//   - GPU host round-trip (H2D copy + kernel + D2H copy -- what you pay for a
//     one-shot call on host data), and
//   - CPU batch (morton/batch.hpp) for reference.
//
// The transfer-bound result is the honest one for one-shot calls; the
// device-only result is what matters when codes live on the GPU across a
// pipeline (the usual reason to be on the GPU at all).

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>

#include "morton/batch.hpp"
#include "morton/morton.hpp"
#include "morton_cuda/morton_cuda.cuh"

using morton::Morton;
using clock_t_ = std::chrono::steady_clock;

int main() {
    using M = Morton<2, 32>;
    using C = M::coord_type;
    using K = M::code_type;
    const std::size_t N = 1u << 25;  // ~33.5M points
    std::mt19937_64 rng(1);
    std::vector<C> x(N), y(N);
    for (std::size_t i = 0; i < N; ++i) { x[i] = C(rng()); y[i] = C(rng()); }
    std::vector<K> out(N);

    std::printf("2D 32-bit encode, N=%zu (%.0f M points)\n", N, N / 1e6);

    // ---- GPU device-only (exclude transfers) ----
    C *dx, *dy; K* dout;
    morton::cuda::check(cudaMalloc(&dx, N * sizeof(C)), "malloc");
    morton::cuda::check(cudaMalloc(&dy, N * sizeof(C)), "malloc");
    morton::cuda::check(cudaMalloc(&dout, N * sizeof(K)), "malloc");
    cudaMemcpy(dx, x.data(), N * sizeof(C), cudaMemcpyHostToDevice);
    cudaMemcpy(dy, y.data(), N * sizeof(C), cudaMemcpyHostToDevice);

    cudaEvent_t e0, e1;
    cudaEventCreate(&e0);
    cudaEventCreate(&e1);
    morton::cuda::encode2_device<32>(dx, dy, dout, N);  // warmup
    cudaDeviceSynchronize();
    const int reps = 20;
    cudaEventRecord(e0);
    for (int r = 0; r < reps; ++r) morton::cuda::encode2_device<32>(dx, dy, dout, N);
    cudaEventRecord(e1);
    cudaEventSynchronize(e1);
    float ms = 0;
    cudaEventElapsedTime(&ms, e0, e1);
    double dev_s = (ms / 1e3) / reps;
    std::printf("  GPU device-only       %8.1f Mops/s  (%.2f ms)\n", N / dev_s / 1e6, dev_s * 1e3);

    cudaFree(dx); cudaFree(dy); cudaFree(dout);

    // ---- GPU host round-trip (alloc+copy+kernel+copy+free) ----
    morton::cuda::encode2_host<32>(x.data(), y.data(), out.data(), N);  // warmup
    auto t0 = clock_t_::now();
    for (int r = 0; r < 5; ++r) morton::cuda::encode2_host<32>(x.data(), y.data(), out.data(), N);
    auto t1 = clock_t_::now();
    double rt_s = std::chrono::duration<double>(t1 - t0).count() / 5;
    std::printf("  GPU host round-trip   %8.1f Mops/s  (%.2f ms, incl. H2D+D2H)\n",
                N / rt_s / 1e6, rt_s * 1e3);

    // ---- CPU batch ----
    morton::batch::encode2<32>(x.data(), y.data(), out.data(), N);  // warmup
    t0 = clock_t_::now();
    for (int r = 0; r < 5; ++r) morton::batch::encode2<32>(x.data(), y.data(), out.data(), N);
    t1 = clock_t_::now();
    double cpu_s = std::chrono::duration<double>(t1 - t0).count() / 5;
    std::printf("  CPU batch (1 core)    %8.1f Mops/s  (%.2f ms)\n", N / cpu_s / 1e6, cpu_s * 1e3);

    return 0;
}
