#include "doctest.h"
#include "morton/iterate.hpp"
#include "morton/morton.hpp"

#include <algorithm>
#include <array>
#include <set>
#include <vector>

using namespace morton;

TEST_CASE("face_neighbors: 2*Dim von Neumann neighbours, exhaustive 2D") {
    using M = Morton<2, 6>;
    for (std::uint32_t x = 1; x < 63; ++x)
        for (std::uint32_t y = 1; y < 63; ++y) {
            M m = M::encode(std::uint8_t(x), std::uint8_t(y));
            auto f = m.face_neighbors();
            REQUIRE(f.size() == 4);
            CHECK(f[0].decode() == std::array<M::coord_type, 2>{std::uint8_t(x - 1), std::uint8_t(y)});
            CHECK(f[1].decode() == std::array<M::coord_type, 2>{std::uint8_t(x + 1), std::uint8_t(y)});
            CHECK(f[2].decode() == std::array<M::coord_type, 2>{std::uint8_t(x), std::uint8_t(y - 1)});
            CHECK(f[3].decode() == std::array<M::coord_type, 2>{std::uint8_t(x), std::uint8_t(y + 1)});
        }
}

TEST_CASE("all_neighbors: Moore neighbourhood sizes and contents") {
    CHECK(Morton<2, 8>::encode(std::uint8_t(5), std::uint8_t(5)).all_neighbors().size() == 8);
    CHECK(Morton<3, 8>::encode(std::uint8_t(5), std::uint8_t(5), std::uint8_t(5))
              .all_neighbors()
              .size() == 26);

    using M = Morton<2, 8>;
    M m = M::encode(std::uint8_t(10), std::uint8_t(20));
    std::set<std::pair<int, int>> got;
    for (auto n : m.all_neighbors()) {
        auto c = n.decode();
        got.insert({c[0], c[1]});
    }
    std::set<std::pair<int, int>> want;
    for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
            if (dx || dy) want.insert({10 + dx, 20 + dy});
    CHECK(got == want);
}

TEST_CASE("ancestor / child_index / child hierarchy") {
    using M = Morton<2, 8>;
    M m = M::encode(std::uint8_t(0b10110101), std::uint8_t(0b01101100));

    // ancestor clears the low level*Dim bits
    for (unsigned level = 0; level <= 8; ++level) {
        M a = m.ancestor(level);
        M::code_type lowmask =
            (level * 2 >= M::code_bits) ? M::field_mask : ((M::code_type(1) << (level * 2)) - 1);
        CHECK((a.code() & lowmask) == 0);
        CHECK((a.code() & ~lowmask) == (m.code() & ~lowmask));
    }

    // child of an ancestor reproduces the original child index
    for (unsigned level = 1; level <= 8; ++level) {
        unsigned idx = m.child_index(level - 1);
        M parent = m.ancestor(level);
        M c = parent.child(level, idx);
        CHECK(c.ancestor(level - 1) == m.ancestor(level - 1));
    }
}

TEST_CASE("range-query helpers find correct codes") {
    using M = Morton<2, 5>;
    using coord = M::coord_type;
    std::array<coord, 2> lo{4, 6}, hi{20, 18};

    std::vector<std::uint32_t> in;
    for (coord y = lo[1]; y <= hi[1]; ++y)
        for (coord x = lo[0]; x <= hi[0]; ++x) in.push_back(M::encode(x, y).code());
    std::sort(in.begin(), in.end());

    for (std::uint32_t z = 0; z < (1u << 10); ++z) {
        M out;
        bool ok = bigmin_in_box<2, 5>(lo, hi, M::from_code(z), out);
        auto up = std::upper_bound(in.begin(), in.end(), z);
        if (up == in.end()) {
            CHECK_FALSE(ok);
        } else {
            REQUIRE(ok);
            CHECK(out.code() == *up);
        }

        ok = litmax_in_box<2, 5>(lo, hi, M::from_code(z), out);
        auto lb = std::lower_bound(in.begin(), in.end(), z);
        if (lb == in.begin()) {
            CHECK_FALSE(ok);
        } else {
            REQUIRE(ok);
            CHECK(out.code() == *(lb - 1));
        }
    }
}
