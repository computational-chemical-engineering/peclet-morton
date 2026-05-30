// Correctness test for the CUDA backend: every GPU result must match the CPU
// Morton library bit-for-bit. Exits non-zero on any mismatch.

#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>

#include "morton/morton.hpp"
#include "morton_cuda/morton_cuda.cuh"

using morton::Morton;

static int failures = 0;
static void expect(bool ok, const char* what) {
    if (!ok) {
        ++failures;
        std::printf("  FAIL: %s\n", what);
    } else {
        std::printf("  ok:   %s\n", what);
    }
}

int main() {
    const std::size_t N = 300000;
    std::mt19937_64 rng(2024);

    // ---- 2D 32-bit encode ----
    std::vector<std::uint32_t> x(N), y(N);
    for (std::size_t i = 0; i < N; ++i) {
        x[i] = std::uint32_t(rng());
        y[i] = std::uint32_t(rng());
    }
    std::vector<std::uint64_t> code(N);
    morton::cuda::encode2_host<32>(x.data(), y.data(), code.data(), N);
    bool ok = true;
    for (std::size_t i = 0; i < N; ++i)
        if (code[i] != Morton<2, 32>::encode(x[i], y[i]).code()) { ok = false; break; }
    expect(ok, "encode2<32> matches CPU");

    // ---- 2D 32-bit decode round-trip ----
    std::vector<std::uint32_t> x2(N), y2(N);
    morton::cuda::decode2_host<32>(code.data(), x2.data(), y2.data(), N);
    ok = (x2 == x && y2 == y);
    expect(ok, "decode2<32> round-trips");

    // ---- 3D 21-bit encode ----
    const std::uint32_t m21 = (1u << 21) - 1;
    std::vector<std::uint32_t> a(N), b(N), c(N);
    for (std::size_t i = 0; i < N; ++i) {
        a[i] = std::uint32_t(rng()) & m21;
        b[i] = std::uint32_t(rng()) & m21;
        c[i] = std::uint32_t(rng()) & m21;
    }
    std::vector<std::uint64_t> code3(N);
    morton::cuda::encode3_host<21>(a.data(), b.data(), c.data(), code3.data(), N);
    ok = true;
    for (std::size_t i = 0; i < N; ++i)
        if (code3[i] != Morton<3, 21>::encode(a[i], b[i], c[i]).code()) { ok = false; break; }
    expect(ok, "encode3<21> matches CPU");

    // ---- per-axis add (positive and negative delta) ----
    std::vector<std::uint64_t> out(N);
    morton::cuda::add_host<2, 32>(code.data(), out.data(), N, /*axis*/ 0, /*delta*/ +5);
    ok = true;
    for (std::size_t i = 0; i < N; ++i) {
        Morton<2, 32> m = Morton<2, 32>::from_code(code[i]);
        m.add(0, 5);
        if (out[i] != m.code()) { ok = false; break; }
    }
    expect(ok, "add<2,32>(axis=0,+5) matches CPU");

    morton::cuda::add_host<2, 32>(code.data(), out.data(), N, /*axis*/ 1, /*delta*/ -3);
    ok = true;
    for (std::size_t i = 0; i < N; ++i) {
        Morton<2, 32> m = Morton<2, 32>::from_code(code[i]);
        m.sub(1, 3);
        if (out[i] != m.code()) { ok = false; break; }
    }
    expect(ok, "add<2,32>(axis=1,-3) matches CPU");

    std::printf(failures ? "\nCUDA TESTS FAILED (%d)\n" : "\nCUDA TESTS PASSED\n", failures);
    return failures ? 1 : 0;
}
