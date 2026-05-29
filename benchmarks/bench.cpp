// Micro-benchmarks for the morton-arithmetic library.
//
// The headline comparison is "walk to the neighbouring cell": doing it with
// the library's O(1) axis arithmetic versus the conventional
// decode -> modify -> re-encode round trip that any encode/decode-only library
// (e.g. libmorton) forces you into.
//
// Build: handled by CMake (-O3). Run: ./morton_bench

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>

#include "morton/iterate.hpp"
#include "morton/morton.hpp"
#include "libmorton/morton.h"

using clk = std::chrono::steady_clock;

template <typename F>
static double time_op(std::size_t iters, F&& f) {
    auto t0 = clk::now();
    f();
    auto t1 = clk::now();
    return std::chrono::duration<double>(t1 - t0).count();
}

static void line(const char* name, std::size_t ops, double secs) {
    printf("  %-46s %8.1f Mops/s  (%6.2f ns/op)\n", name, ops / secs / 1e6,
           secs / ops * 1e9);
}

int main() {
    using M = morton::Morton<2, 32>;
    std::mt19937_64 rng(42);

    const std::size_t N = 1u << 24;  // ~16.7M
    std::vector<std::uint32_t> xs(N), ys(N);
    for (std::size_t i = 0; i < N; ++i) {
        xs[i] = std::uint32_t(rng());
        ys[i] = std::uint32_t(rng());
    }
    volatile std::uint64_t sink = 0;

    printf("== 2D 32-bit encode/decode (N=%zu) ==\n", N);

    double t = time_op(N, [&] {
        std::uint64_t acc = 0;
        for (std::size_t i = 0; i < N; ++i) acc ^= M::encode(xs[i], ys[i]).code();
        sink = acc;
    });
    line("encode (this library, PDEP)", N, t);

    t = time_op(N, [&] {
        std::uint64_t acc = 0;
        for (std::size_t i = 0; i < N; ++i) acc ^= libmorton::morton2D_64_encode(xs[i], ys[i]);
        sink = acc;
    });
    line("encode (libmorton)", N, t);

    std::vector<std::uint64_t> codes(N);
    for (std::size_t i = 0; i < N; ++i) codes[i] = M::encode(xs[i], ys[i]).code();

    t = time_op(N, [&] {
        std::uint64_t acc = 0;
        for (std::size_t i = 0; i < N; ++i) {
            auto d = M::from_code(codes[i]).decode();
            acc ^= d[0] ^ d[1];
        }
        sink = acc;
    });
    line("decode (this library, PEXT)", N, t);

    t = time_op(N, [&] {
        std::uint64_t acc = 0;
        for (std::size_t i = 0; i < N; ++i) {
            std::uint_fast32_t x, y;
            libmorton::morton2D_64_decode(codes[i], x, y);
            acc ^= x ^ y;
        }
        sink = acc;
    });
    line("decode (libmorton)", N, t);

    printf("\n== Neighbour step: +1 along axis 0 (N=%zu) ==\n", N);

    // The library way: O(1) arithmetic, code stays interleaved.
    t = time_op(N, [&] {
        std::uint64_t acc = 0;
        for (std::size_t i = 0; i < N; ++i) {
            M m = M::from_code(codes[i]);
            m.inc(0);
            acc ^= m.code();
        }
        sink = acc;
    });
    line("neighbour via arithmetic (this library)", N, t);

    // The conventional way with an encode/decode-only library: round trip.
    t = time_op(N, [&] {
        std::uint64_t acc = 0;
        for (std::size_t i = 0; i < N; ++i) {
            std::uint_fast32_t x, y;
            libmorton::morton2D_64_decode(codes[i], x, y);
            acc ^= libmorton::morton2D_64_encode(std::uint_fast32_t(x + 1), y);
        }
        sink = acc;
    });
    line("neighbour via decode+encode (libmorton)", N, t);

    printf("\n== 4-neighbour stencil on scattered codes (N=%zu) ==\n", N);
    // For each code, combine its four axis-neighbours (-x,+x,-y,+y). This is
    // the access pattern of a stencil/finite-difference sweep over a Z-ordered
    // grid, where you hold a code and need its neighbours.
    t = time_op(N, [&] {
        std::uint64_t acc = 0;
        for (std::size_t i = 0; i < N; ++i) {
            M c = M::from_code(codes[i]);
            M xm = c; xm.dec(0);
            M xp = c; xp.inc(0);
            M ym = c; ym.dec(1);
            M yp = c; yp.inc(1);
            acc ^= xm.code() ^ xp.code() ^ ym.code() ^ yp.code();
        }
        sink = acc;
    });
    line("stencil via arithmetic (this library)", N * 4, t);

    t = time_op(N, [&] {
        std::uint64_t acc = 0;
        for (std::size_t i = 0; i < N; ++i) {
            std::uint_fast32_t x, y;
            libmorton::morton2D_64_decode(codes[i], x, y);
            acc ^= libmorton::morton2D_64_encode(std::uint_fast32_t(x - 1), y);
            acc ^= libmorton::morton2D_64_encode(std::uint_fast32_t(x + 1), y);
            acc ^= libmorton::morton2D_64_encode(x, std::uint_fast32_t(y - 1));
            acc ^= libmorton::morton2D_64_encode(x, std::uint_fast32_t(y + 1));
        }
        sink = acc;
    });
    line("stencil via decode+encode (libmorton)", N * 4, t);

    printf("\n== Dense region sweep: 4096x4096 cells ==\n");
    const std::uint32_t W = 4096;
    const std::size_t cells = std::size_t(W) * W;

    // Row-major arithmetic walk (never re-encodes).
    t = time_op(cells, [&] {
        std::uint64_t acc = 0;
        morton::for_each_in_box<2, 32>({1000, 2000}, {1000 + W - 1, 2000 + W - 1},
                                       [&](M m) { acc ^= m.code(); });
        sink = acc;
    });
    line("sweep via for_each_in_box (arithmetic)", cells, t);

    // Naive: encode each (x,y) from scratch.
    t = time_op(cells, [&] {
        std::uint64_t acc = 0;
        for (std::uint32_t y = 2000; y < 2000 + W; ++y)
            for (std::uint32_t x = 1000; x < 1000 + W; ++x)
                acc ^= M::encode(x, y).code();
        sink = acc;
    });
    line("sweep via per-cell encode (PDEP)", cells, t);

    (void)sink;
    return 0;
}
