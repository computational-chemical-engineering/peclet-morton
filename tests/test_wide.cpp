#include "doctest.h"
#include "morton/morton.hpp"

#include <array>
#include <cstdint>
#include <random>

using namespace morton;

#if defined(MORTON_HAS_INT128)

// Reference bit-by-bit encoder over 128-bit codes.
template <unsigned Dim, unsigned Bits, typename Code>
static Code ref_encode(const std::array<std::uint64_t, Dim>& c) {
    Code m = 0;
    for (unsigned b = 0; b < Bits; ++b)
        for (unsigned d = 0; d < Dim; ++d)
            m |= (Code((c[d] >> b) & 1ull)) << (b * Dim + d);
    return m;
}

TEST_CASE("2D 64-bit (128-bit code) round-trips and arithmetic") {
    using M = Morton<2, 64>;
    static_assert(M::code_bits == 128, "");
    static_assert(sizeof(M::code_type) == 16, "");

    std::mt19937_64 rng(1);
    for (int i = 0; i < 20000; ++i) {
        std::array<std::uint64_t, 2> c{rng(), rng()};
        M m = M::encode(c[0], c[1]);
        CHECK(m.code() == ref_encode<2, 64, M::code_type>(c));
        auto d = m.decode();
        CHECK(d[0] == c[0]);
        CHECK(d[1] == c[1]);

        unsigned ax = rng() & 1;
        std::uint64_t k = rng();
        M a = m; a.add(ax, k);
        auto da = a.decode();
        CHECK(da[ax] == std::uint64_t(c[ax] + k));   // wraps mod 2^64
        CHECK(da[ax ^ 1] == c[ax ^ 1]);
    }
}

TEST_CASE("3D 32-bit (96-bit code) round-trips and arithmetic") {
    using M = Morton<3, 32>;
    static_assert(M::code_bits == 96, "");

    std::mt19937_64 rng(2);
    for (int i = 0; i < 20000; ++i) {
        std::array<std::uint64_t, 3> c{std::uint32_t(rng()), std::uint32_t(rng()),
                                       std::uint32_t(rng())};
        M m = M::encode(std::uint32_t(c[0]), std::uint32_t(c[1]), std::uint32_t(c[2]));
        CHECK(m.code() == ref_encode<3, 32, M::code_type>(c));
        auto d = m.decode();
        CHECK(d[0] == c[0]);
        CHECK(d[1] == c[1]);
        CHECK(d[2] == c[2]);

        M n = m.neighbor(2, +1);
        CHECK(n.get(2) == std::uint32_t(c[2] + 1));
        CHECK(n.get(0) == c[0]);
        CHECK(n.get(1) == c[1]);
    }
}

TEST_CASE("128-bit field mask is full width") {
    CHECK(Morton<2, 64>::field_mask == Morton<2, 64>::code_type(~Morton<2, 64>::code_type(0)));
}

#else
TEST_CASE("128-bit codes unavailable on this compiler") {
    WARN_MESSAGE(false, "no __int128: wide-code tests skipped");
}
#endif

// ---- codes wider than 128 bits use the wide_uint<W> software backend -------

// Reference encoder working in the library's own wide code type.
template <unsigned Dim, unsigned Bits, typename Code>
static Code wide_ref_encode(const std::array<std::uint64_t, Dim>& c) {
    Code m = 0;
    for (unsigned b = 0; b < Bits; ++b)
        for (unsigned d = 0; d < Dim; ++d) {
            Code bit = Code((c[d] >> b) & 1ull);
            m = m | (bit << (b * Dim + d));
        }
    return m;
}

TEST_CASE("3D 64-bit (192-bit code) round-trips and arithmetic") {
    using M = Morton<3, 64>;
    static_assert(M::code_bits == 192, "");
    static_assert(sizeof(M::code_type) == 24, "wide_uint<3>");

    std::mt19937_64 rng(101);
    for (int i = 0; i < 5000; ++i) {
        std::array<std::uint64_t, 3> c{rng(), rng(), rng()};
        M m = M::encode(c[0], c[1], c[2]);
        CHECK(m.code() == wide_ref_encode<3, 64, M::code_type>(c));
        auto d = m.decode();
        CHECK(d[0] == c[0]);
        CHECK(d[1] == c[1]);
        CHECK(d[2] == c[2]);

        unsigned ax = rng() % 3;
        std::uint64_t k = rng();
        M a = m; a.add(ax, k);
        auto da = a.decode();
        CHECK(da[ax] == std::uint64_t(c[ax] + k));  // wraps mod 2^64
        CHECK(da[(ax + 1) % 3] == c[(ax + 1) % 3]);
        CHECK(da[(ax + 2) % 3] == c[(ax + 2) % 3]);
    }
}

TEST_CASE("2D 128-bit (256-bit code) round-trips, neighbours, Z-order") {
    using M = Morton<2, 128>;
    static_assert(M::code_bits == 256, "");
    static_assert(sizeof(M::code_type) == 32, "wide_uint<4>");
    static_assert(M::field_mask == M::code_type(~M::code_type(0)), "");

#if defined(MORTON_HAS_INT128)
    std::mt19937_64 rng(202);
    for (int i = 0; i < 5000; ++i) {
        unsigned __int128 x = ((unsigned __int128)(rng()) << 64) | rng();
        unsigned __int128 y = ((unsigned __int128)(rng()) << 64) | rng();
        M m = M::encode(x, y);
        auto d = m.decode();
        CHECK(d[0] == x);
        CHECK(d[1] == y);

        M n = m.neighbor(0, +1);
        CHECK(n.get(0) == (unsigned __int128)(x + 1));
        CHECK(n.get(1) == y);
    }
#endif
}

TEST_CASE("wide-code exhaustive small-value arithmetic matches built-in path") {
    // Compare a 192-bit code restricted to small coordinates against the same
    // logical operation: an independent recompute.
    using M = Morton<3, 64>;
    for (std::uint64_t x = 0; x < 40; ++x)
        for (std::uint64_t y = 0; y < 40; ++y) {
            M m = M::encode(x, y, std::uint64_t(7));
            m.inc(0);
            m.dec(1);
            auto c = m.decode();
            REQUIRE(c[0] == x + 1);
            REQUIRE(c[1] == (y == 0 ? std::uint64_t(~0ull) : y - 1));
            REQUIRE(c[2] == 7);
        }
}
