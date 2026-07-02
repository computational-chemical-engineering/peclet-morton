#include <array>
#include <cstdint>
#include <random>

#include "doctest.h"
#include "libmorton/morton.h"
#include "morton/morton.hpp"

using namespace morton;

// Independent, obviously-correct reference encoder/decoder (bit by bit).
template <unsigned Dim, unsigned Bits>
static std::uint64_t ref_encode(const std::array<std::uint64_t, Dim>& c) {
  std::uint64_t m = 0;
  for (unsigned b = 0; b < Bits; ++b)
    for (unsigned d = 0; d < Dim; ++d)
      m |= ((c[d] >> b) & 1ull) << (b * Dim + d);
  return m;
}

TEST_CASE("encode/decode round-trips against the bit-by-bit reference") {
  std::mt19937_64 rng(12345);

  SUBCASE("2D 32-bit") {
    using M = Morton<2, 32>;
    for (int i = 0; i < 20000; ++i) {
      std::array<std::uint64_t, 2> c{rng() & 0xffffffffu, rng() & 0xffffffffu};
      M m = M::encode(std::uint32_t(c[0]), std::uint32_t(c[1]));
      CHECK(m.code() == ref_encode<2, 32>(c));
      auto d = m.decode();
      CHECK(d[0] == c[0]);
      CHECK(d[1] == c[1]);
    }
  }
  SUBCASE("3D 21-bit") {
    using M = Morton<3, 21>;
    const std::uint64_t mask = (1ull << 21) - 1;
    for (int i = 0; i < 20000; ++i) {
      std::array<std::uint64_t, 3> c{rng() & mask, rng() & mask, rng() & mask};
      M m = M::encode(std::uint32_t(c[0]), std::uint32_t(c[1]), std::uint32_t(c[2]));
      CHECK(m.code() == ref_encode<3, 21>(c));
      auto d = m.decode();
      CHECK(d[0] == c[0]);
      CHECK(d[1] == c[1]);
      CHECK(d[2] == c[2]);
    }
  }
  SUBCASE("3D 16-bit and 2D 16-bit") {
    using M3 = Morton<3, 16>;
    using M2 = Morton<2, 16>;
    for (int i = 0; i < 20000; ++i) {
      std::array<std::uint64_t, 3> c{rng() & 0xffff, rng() & 0xffff, rng() & 0xffff};
      CHECK(M3::encode(std::uint16_t(c[0]), std::uint16_t(c[1]), std::uint16_t(c[2])).code() ==
            ref_encode<3, 16>(c));
      std::array<std::uint64_t, 2> c2{c[0], c[1]};
      CHECK(M2::encode(std::uint16_t(c2[0]), std::uint16_t(c2[1])).code() == ref_encode<2, 16>(c2));
    }
  }
}

TEST_CASE("agreement with libmorton (the de-facto reference library)") {
  std::mt19937_64 rng(999);
  SUBCASE("2D 32-bit") {
    for (int i = 0; i < 20000; ++i) {
      std::uint32_t x = std::uint32_t(rng()), y = std::uint32_t(rng());
      CHECK(Morton<2, 32>::encode(x, y).code() == libmorton::morton2D_64_encode(x, y));
    }
  }
  SUBCASE("3D 21-bit") {
    const std::uint32_t mask = (1u << 21) - 1;
    for (int i = 0; i < 20000; ++i) {
      std::uint32_t x = std::uint32_t(rng()) & mask;
      std::uint32_t y = std::uint32_t(rng()) & mask;
      std::uint32_t z = std::uint32_t(rng()) & mask;
      CHECK(Morton<3, 21>::encode(x, y, z).code() == libmorton::morton3D_64_encode(x, y, z));
    }
  }
}

TEST_CASE("exhaustive round-trip for small grids") {
  using M = Morton<2, 6>;  // 64x64
  for (std::uint32_t x = 0; x < 64; ++x)
    for (std::uint32_t y = 0; y < 64; ++y) {
      auto d = M::encode(std::uint8_t(x), std::uint8_t(y)).decode();
      REQUIRE(d[0] == x);
      REQUIRE(d[1] == y);
    }
}

TEST_CASE("constants and masks") {
  CHECK(Morton<2, 32>::field_mask == 0xffffffffffffffffull);
  CHECK(Morton<3, 21>::field_mask == ((1ull << 63) - 1));
  CHECK(Morton<2, 8>::coord_max == 255);
  // axis masks partition the field with no overlap
  using M = Morton<3, 10>;
  auto m = M::axis_mask(0) | M::axis_mask(1) | M::axis_mask(2);
  CHECK(m == M::field_mask);
  CHECK((M::axis_mask(0) & M::axis_mask(1)) == 0);
  CHECK((M::axis_mask(1) & M::axis_mask(2)) == 0);
}
