// Correctness test for the portable Kokkos backend: every device result must
// match the scalar Morton library bit-for-bit, on whichever backend Kokkos was
// built with (CUDA / HIP / OpenMP / Serial). Exits non-zero on any mismatch.
//
// This replaces cuda/tests/test_cuda.cu and extends coverage to all four
// binding-supported layouts: (2,32) (2,16) (3,21) (3,16).

#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>

#include "morton/kokkos.hpp"
#include "morton/morton.hpp"

using morton::Morton;
namespace mk = morton::kokkos;

static int failures = 0;
static void expect(bool ok, const char* what) {
    if (!ok) {
        ++failures;
        std::printf("  FAIL: %s\n", what);
    } else {
        std::printf("  ok:   %s\n", what);
    }
}

// Encode/decode round-trip + scalar cross-check for a 2D layout.
template <unsigned Bits>
static void check2d(std::size_t N, std::mt19937_64& rng) {
    using M = Morton<2, Bits>;
    using C = typename M::coord_type;
    using K = typename M::code_type;
    const C mask = M::coord_max;
    std::vector<C> x(N), y(N);
    for (std::size_t i = 0; i < N; ++i) {
        x[i] = C(rng()) & mask;
        y[i] = C(rng()) & mask;
    }
    std::vector<K> code(N);
    mk::encode2_host<Bits>(x.data(), y.data(), code.data(), N);
    bool ok = true;
    for (std::size_t i = 0; i < N; ++i)
        if (code[i] != M::encode(x[i], y[i]).code()) { ok = false; break; }
    expect(ok, "encode2 matches scalar");

    std::vector<C> x2(N), y2(N);
    mk::decode2_host<Bits>(code.data(), x2.data(), y2.data(), N);
    expect(x2 == x && y2 == y, "decode2 round-trips");
}

// Encode/decode round-trip + scalar cross-check for a 3D layout.
template <unsigned Bits>
static void check3d(std::size_t N, std::mt19937_64& rng) {
    using M = Morton<3, Bits>;
    using C = typename M::coord_type;
    using K = typename M::code_type;
    const C mask = M::coord_max;
    std::vector<C> a(N), b(N), c(N);
    for (std::size_t i = 0; i < N; ++i) {
        a[i] = C(rng()) & mask;
        b[i] = C(rng()) & mask;
        c[i] = C(rng()) & mask;
    }
    std::vector<K> code(N);
    mk::encode3_host<Bits>(a.data(), b.data(), c.data(), code.data(), N);
    bool ok = true;
    for (std::size_t i = 0; i < N; ++i)
        if (code[i] != M::encode(a[i], b[i], c[i]).code()) { ok = false; break; }
    expect(ok, "encode3 matches scalar");

    std::vector<C> a2(N), b2(N), c2(N);
    mk::decode3_host<Bits>(code.data(), a2.data(), b2.data(), c2.data(), N);
    expect(a2 == a && b2 == b && c2 == c, "decode3 round-trips");
}

int main(int argc, char** argv) {
    Kokkos::ScopeGuard guard(argc, argv);
    const std::size_t N = 300000;
    std::mt19937_64 rng(2024);

    std::printf("(2,32):\n");
    check2d<32>(N, rng);
    std::printf("(2,16):\n");
    check2d<16>(N, rng);
    std::printf("(3,21):\n");
    check3d<21>(N, rng);
    std::printf("(3,16):\n");
    check3d<16>(N, rng);

    // ---- per-axis add (positive and negative delta), wraps mod 2^Bits ----
    {
        using M = Morton<2, 32>;
        std::vector<std::uint32_t> x(N), y(N);
        for (std::size_t i = 0; i < N; ++i) {
            x[i] = std::uint32_t(rng());
            y[i] = std::uint32_t(rng());
        }
        std::vector<std::uint64_t> code(N), out(N);
        mk::encode2_host<32>(x.data(), y.data(), code.data(), N);

        mk::add_host<2, 32>(code.data(), out.data(), N, /*axis*/ 0, /*delta*/ +5);
        bool ok = true;
        for (std::size_t i = 0; i < N; ++i) {
            M m = M::from_code(code[i]);
            m.add(0, 5);
            if (out[i] != m.code()) { ok = false; break; }
        }
        expect(ok, "add<2,32>(axis=0,+5) matches scalar");

        mk::add_host<2, 32>(code.data(), out.data(), N, /*axis*/ 1, /*delta*/ -3);
        ok = true;
        for (std::size_t i = 0; i < N; ++i) {
            M m = M::from_code(code[i]);
            m.sub(1, 3);
            if (out[i] != m.code()) { ok = false; break; }
        }
        expect(ok, "add<2,32>(axis=1,-3) matches scalar");
    }

    std::printf(failures ? "\nKOKKOS TESTS FAILED (%d)\n" : "\nKOKKOS TESTS PASSED\n", failures);
    return failures ? 1 : 0;
}
