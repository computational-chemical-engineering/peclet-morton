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
