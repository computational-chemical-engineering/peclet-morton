// Batch / SIMD benchmark. The per-axis arithmetic auto-vectorises; this shows
// the throughput when it is cache-resident (vectorisation visible) vs. when the
// arrays exceed cache (memory-bound, the honest common case).

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>

#include "morton/batch.hpp"

using clk = std::chrono::steady_clock;
using morton::Morton;

template <typename F>
static double best_secs(int reps, F&& f) {
  double b = 1e30;
  for (int r = 0; r < reps; ++r) {
    auto t0 = clk::now();
    f();
    auto t1 = clk::now();
    b = std::min(b, std::chrono::duration<double>(t1 - t0).count());
  }
  return b;
}

int main() {
  using M = Morton<2, 32>;
  using code = M::code_type;
  std::mt19937_64 rng(7);
  volatile code sink = 0;

#if MORTON_X86
  // Report which path morton::batch will take on this CPU. The AVX-512 batch
  // kernels only show their advantage on hardware that has AVX-512F; on other
  // CPUs (and under emulation) the auto-vectorised scalar path runs instead.
  printf("CPU dispatch: BMI2=%d AVX2=%d AVX-512F=%d  (batch encode uses %s)\n",
         int(morton::detail::cpu_has_bmi2()), int(morton::detail::cpu_has_avx2()),
         int(morton::detail::cpu_has_avx512f()),
         morton::detail::cpu_has_avx512f() ? "AVX-512" : "scalar/AVX2");
#endif

  for (std::size_t N : {std::size_t(1u << 12), std::size_t(1u << 24)}) {
    std::vector<code> in(N), out(N);
    for (auto& v : in)
      v = M::encode(std::uint32_t(rng()), std::uint32_t(rng())).code();
    const char* tag = (N <= (1u << 16)) ? "cache-resident" : "out-of-cache";

    double t = best_secs(50, [&] { morton::batch::add<2, 32>(in.data(), out.data(), N, 0, 1); });
    sink = out[N - 1];
    printf("N=%-9zu (%-14s) batch::add(axis0,+1)  %8.1f Mops/s\n", N, tag, N / t / 1e6);

    t = best_secs(50, [&] {
      for (std::size_t i = 0; i < N; ++i) {
        M m = M::from_code(in[i]);
        m.inc(0);
        out[i] = m.code();
      }
    });
    sink = out[N - 1];
    printf("N=%-9zu (%-14s) scalar per-object inc %8.1f Mops/s\n", N, tag, N / t / 1e6);
  }

  // Explicit AVX-512 batch encode vs. per-element encode. The bulk encode is
  // the cache-resident transform the AVX-512 kernel targets; on AVX-512
  // hardware morton::batch::encode2 runs the vector kernel.
  for (std::size_t N : {std::size_t(1u << 12), std::size_t(1u << 24)}) {
    std::vector<std::uint32_t> x(N), y(N);
    std::vector<code> out(N);
    for (std::size_t i = 0; i < N; ++i) {
      x[i] = std::uint32_t(rng());
      y[i] = std::uint32_t(rng());
    }
    const char* tag = (N <= (1u << 16)) ? "cache-resident" : "out-of-cache";

    double t =
        best_secs(50, [&] { morton::batch::encode2<32>(x.data(), y.data(), out.data(), N); });
    sink = out[N - 1];
    printf("N=%-9zu (%-14s) batch::encode2        %8.1f Mops/s\n", N, tag, N / t / 1e6);

    t = best_secs(50, [&] {
      for (std::size_t i = 0; i < N; ++i)
        out[i] = M::encode(x[i], y[i]).code();
    });
    sink = out[N - 1];
    printf("N=%-9zu (%-14s) scalar per-elem encode%8.1f Mops/s\n", N, tag, N / t / 1e6);
  }
  (void)sink;
  return 0;
}
