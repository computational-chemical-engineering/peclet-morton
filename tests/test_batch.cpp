#include <cstdint>
#include <random>
#include <vector>

#include "doctest.h"
#include "morton/batch.hpp"
#include "morton/morton.hpp"

using namespace morton;

TEST_CASE("batch::add/sub/step match the scalar object API") {
  using M = Morton<2, 32>;
  using code = M::code_type;
  std::mt19937_64 rng(11);
  const std::size_t n = 5000;
  std::vector<code> in(n), out(n);
  for (auto& v : in)
    v = M::encode(std::uint32_t(rng()), std::uint32_t(rng())).code();

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
  for (std::size_t i = 0; i < n; ++i)
    CHECK(out[i] == M::encode(x[i], y[i], z[i]).code());
}

// These exercise the AVX-512 dispatch in batch.hpp. On a CPU with AVX-512F the
// vector kernels run; otherwise the scalar fallback. Either way the result must
// equal the per-element core path. Sizes are deliberately NOT multiples of 8 so
// the AVX-512 scalar tail (the i<n remainder loop) is covered. Run the whole
// suite under `sde64 -skx -- ./morton_tests` to validate the AVX-512 path on a
// machine that lacks AVX-512.

TEST_CASE("batch::encode2/decode2 (2,32) round-trip matches scalar, incl. tail") {
  using M = Morton<2, 32>;
  std::mt19937_64 rng(7);
  for (std::size_t n : {std::size_t(0), std::size_t(1), std::size_t(7), std::size_t(8),
                        std::size_t(9), std::size_t(1003)}) {
    std::vector<std::uint32_t> x(n), y(n);
    for (std::size_t i = 0; i < n; ++i) {
      x[i] = std::uint32_t(rng());
      y[i] = std::uint32_t(rng());
    }
    std::vector<M::code_type> code(n);
    batch::encode2<32>(x.data(), y.data(), code.data(), n);
    for (std::size_t i = 0; i < n; ++i)
      REQUIRE(code[i] == M::encode(x[i], y[i]).code());

    std::vector<std::uint32_t> dx(n), dy(n);
    batch::decode2<32>(code.data(), dx.data(), dy.data(), n);
    for (std::size_t i = 0; i < n; ++i) {
      REQUIRE(dx[i] == x[i]);
      REQUIRE(dy[i] == y[i]);
    }
  }
}

TEST_CASE("batch::encode3/decode3 (3,21) round-trip matches scalar, incl. tail") {
  using M = Morton<3, 21>;
  std::mt19937_64 rng(13);
  const std::uint32_t mask = (1u << 21) - 1;
  for (std::size_t n :
       {std::size_t(0), std::size_t(5), std::size_t(8), std::size_t(11), std::size_t(2049)}) {
    std::vector<std::uint32_t> x(n), y(n), z(n);
    for (std::size_t i = 0; i < n; ++i) {
      x[i] = std::uint32_t(rng()) & mask;
      y[i] = std::uint32_t(rng()) & mask;
      z[i] = std::uint32_t(rng()) & mask;
    }
    std::vector<M::code_type> code(n);
    batch::encode3<21>(x.data(), y.data(), z.data(), code.data(), n);
    for (std::size_t i = 0; i < n; ++i)
      REQUIRE(code[i] == M::encode(x[i], y[i], z[i]).code());

    std::vector<std::uint32_t> dx(n), dy(n), dz(n);
    batch::decode3<21>(code.data(), dx.data(), dy.data(), dz.data(), n);
    for (std::size_t i = 0; i < n; ++i) {
      REQUIRE(dx[i] == x[i]);
      REQUIRE(dy[i] == y[i]);
      REQUIRE(dz[i] == z[i]);
    }
  }
}

TEST_CASE("batch::add/sub on 64-bit codes match scalar across the AVX-512 path, incl. tail") {
  using M = Morton<3, 21>;  // 63-bit code -> uint64, hits the AVX-512 add64/sub64
  using code = M::code_type;
  std::mt19937_64 rng(29);
  const std::uint32_t mask = (1u << 21) - 1;
  for (std::size_t n : {std::size_t(0), std::size_t(3), std::size_t(8), std::size_t(1001)}) {
    std::vector<code> in(n), out(n);
    for (auto& v : in)
      v = M::encode(std::uint32_t(rng()) & mask, std::uint32_t(rng()) & mask,
                    std::uint32_t(rng()) & mask)
              .code();
    for (unsigned axis = 0; axis < 3; ++axis) {
      std::uint32_t k = std::uint32_t(rng()) & mask;
      batch::add<3, 21>(in.data(), out.data(), n, axis, k);
      for (std::size_t i = 0; i < n; ++i) {
        M m = M::from_code(in[i]);
        m.add(axis, k);
        REQUIRE(out[i] == m.code());
      }
      batch::sub<3, 21>(in.data(), out.data(), n, axis, k);
      for (std::size_t i = 0; i < n; ++i) {
        M m = M::from_code(in[i]);
        m.sub(axis, k);
        REQUIRE(out[i] == m.code());
      }
    }
  }
}
