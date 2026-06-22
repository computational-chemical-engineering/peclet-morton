// Throughput benchmark for the portable Kokkos backend (replaces the retired
// raw-CUDA bench). Runs on whichever backend Kokkos was built with.
//
// Reports three numbers for 2D 32-bit encode:
//   - device-only      (data already resident in device Views; pure kernel),
//   - host round-trip  (H2D copy + kernel + D2H copy -- what a one-shot call on
//     host data costs), and
//   - CPU batch (morton/batch.hpp) for reference.
//
// The transfer-bound result is the honest one for one-shot calls; the
// device-only result is what matters when codes live on the device across a
// pipeline (the usual reason to be on the GPU at all). On the OpenMP/Serial
// backends the "device" is the host, so the two device numbers converge.

#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>

#include "morton/batch.hpp"
#include "morton/kokkos.hpp"
#include "morton/morton.hpp"

using morton::Morton;
namespace mk = morton::kokkos;

int main(int argc, char** argv) {
    Kokkos::ScopeGuard guard(argc, argv);
    using M = Morton<2, 32>;
    using C = M::coord_type;
    using K = M::code_type;
    const std::size_t N = 1u << 25;  // ~33.5M points
    std::mt19937_64 rng(1);
    std::vector<C> x(N), y(N);
    for (std::size_t i = 0; i < N; ++i) { x[i] = C(rng()); y[i] = C(rng()); }
    std::vector<K> out(N);

    std::printf("2D 32-bit encode, N=%zu (%.0f M points)  [Kokkos default exec space]\n",
                N, N / 1e6);

    // ---- device-only (exclude transfers): data resident in device Views ----
    Kokkos::View<C*> dx("x", N), dy("y", N);
    Kokkos::View<K*> dout("out", N);
    using HU = Kokkos::View<const C*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    Kokkos::deep_copy(dx, HU(x.data(), N));
    Kokkos::deep_copy(dy, HU(y.data(), N));

    mk::encode2<32>(dx, dy, dout);  // warmup
    Kokkos::fence();
    const int reps = 20;
    Kokkos::Timer timer;
    for (int r = 0; r < reps; ++r) mk::encode2<32>(dx, dy, dout);
    Kokkos::fence();
    double dev_s = timer.seconds() / reps;
    std::printf("  device-only           %8.1f Mops/s  (%.2f ms)\n", N / dev_s / 1e6, dev_s * 1e3);

    // ---- host round-trip (stage + kernel + copy back) ----
    mk::encode2_host<32>(x.data(), y.data(), out.data(), N);  // warmup
    timer.reset();
    for (int r = 0; r < 5; ++r) mk::encode2_host<32>(x.data(), y.data(), out.data(), N);
    double rt_s = timer.seconds() / 5;
    std::printf("  host round-trip       %8.1f Mops/s  (%.2f ms, incl. H<->D)\n",
                N / rt_s / 1e6, rt_s * 1e3);

    // ---- CPU batch ----
    morton::batch::encode2<32>(x.data(), y.data(), out.data(), N);  // warmup
    timer.reset();
    for (int r = 0; r < 5; ++r) morton::batch::encode2<32>(x.data(), y.data(), out.data(), N);
    double cpu_s = timer.seconds() / 5;
    std::printf("  CPU batch (1 core)    %8.1f Mops/s  (%.2f ms)\n", N / cpu_s / 1e6, cpu_s * 1e3);

    return 0;
}
