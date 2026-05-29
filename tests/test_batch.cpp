#include "doctest.h"
#include "morton/batch.hpp"
#include "morton/morton.hpp"

#include <cstdint>
#include <random>
#include <vector>

using namespace morton;

TEST_CASE("batch::add/sub/step match the scalar object API") {
    using M = Morton<2, 32>;
    using code = M::code_type;
    std::mt19937_64 rng(11);
    const std::size_t n = 5000;
    std::vector<code> in(n), out(n);
    for (auto& v : in) v = M::encode(std::uint32_t(rng()), std::uint32_t(rng())).code();

    for (unsigned axis = 0; axis < 2; ++axis) {
        std::uint32_t k = std::uint32_t(rng());
        batch::add<2, 32>(in.data(), out.data(), n, axis, k);
        for (std::size_t i = 0; i < n; ++i) {
            M m = M::from_code(in[i]);
            m.add(axis, k);
            REQUIRE(out[i] == m.code());
        }
        batch::sub<2, 32>(in.data(), out.data(), n, axis, k);
        for (std::size_t i = 0; i < n; ++i) {
            M m = M::from_code(in[i]);
            m.sub(axis, k);
            REQUIRE(out[i] == m.code());
        }
        batch::step<2, 32>(in.data(), out.data(), n, axis, +1);
        for (std::size_t i = 0; i < n; ++i) {
            M m = M::from_code(in[i]);
            m.inc(axis);
            REQUIRE(out[i] == m.code());
        }
    }
}

TEST_CASE("batch::add works in place (in == out)") {
    using M = Morton<3, 21>;
    using code = M::code_type;
    std::vector<code> a{M::encode(1u, 2u, 3u).code(), M::encode(10u, 20u, 30u).code()};
    std::vector<code> b = a;
    batch::add<3, 21>(a.data(), a.data(), a.size(), 1, 100);  // in place
    for (std::size_t i = 0; i < b.size(); ++i) {
        M m = M::from_code(b[i]);
        m.add(1, 100);
        CHECK(a[i] == m.code());
    }
}

TEST_CASE("batch encode matches scalar encode") {
    using M = Morton<3, 21>;
    std::mt19937_64 rng(3);
    const std::size_t n = 4000;
    std::vector<std::uint32_t> x(n), y(n), z(n);
    const std::uint32_t mask = (1u << 21) - 1;
    for (std::size_t i = 0; i < n; ++i) {
        x[i] = std::uint32_t(rng()) & mask;
        y[i] = std::uint32_t(rng()) & mask;
        z[i] = std::uint32_t(rng()) & mask;
    }
    std::vector<M::code_type> out(n);
    batch::encode3<21>(x.data(), y.data(), z.data(), out.data(), n);
    for (std::size_t i = 0; i < n; ++i) CHECK(out[i] == M::encode(x[i], y[i], z[i]).code());
}
