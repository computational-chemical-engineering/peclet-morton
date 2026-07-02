#include <array>
#include <cstdint>
#include <random>

#include "doctest.h"
#include "morton/morton.hpp"

using namespace morton;

TEST_CASE("exhaustive 2D 8-bit: arithmetic matches decode/modify/encode") {
  using M = Morton<2, 8>;
  for (std::uint32_t x = 0; x < 256; ++x) {
    for (std::uint32_t y = 0; y < 256; ++y) {
      M base = M::encode(std::uint8_t(x), std::uint8_t(y));

      // inc / dec on each axis, with wrap
      {
        M m = base;
        m.inc(0);
        auto c = m.decode();
        CHECK(c[0] == ((x + 1) & 0xff));
        CHECK(c[1] == y);
      }
      {
        M m = base;
        m.dec(1);
        auto c = m.decode();
        CHECK(c[0] == x);
        CHECK(c[1] == ((y - 1) & 0xff));
      }
      // neighbour helpers
      CHECK(base.neighbor(0, +1) == [&] {
        M m = base;
        m.inc(0);
        return m;
      }());
      CHECK(base.neighbor(1, -1) == [&] {
        M m = base;
        m.dec(1);
        return m;
      }());

      // set replaces one axis, leaves the other
      {
        M m = base;
        m.set(0, 42);
        auto c = m.decode();
        CHECK(c[0] == 42);
        CHECK(c[1] == y);
      }
    }
  }
}

TEST_CASE("add/sub with arbitrary k wrap mod 2^Bits") {
  using M = Morton<3, 8>;
  std::mt19937_64 rng(7);
  for (int i = 0; i < 50000; ++i) {
    std::uint8_t x = rng(), y = rng(), z = rng(), k = rng();
    unsigned d = rng() % 3;
    M m = M::encode(x, y, z);
    m.add(d, k);
    auto c = m.decode();
    std::array<unsigned, 3> exp{x, y, z};
    exp[d] = (exp[d] + k) & 0xff;
    CHECK(c[0] == exp[0]);
    CHECK(c[1] == exp[1]);
    CHECK(c[2] == exp[2]);

    M m2 = M::encode(x, y, z);
    m2.sub(d, k);
    auto c2 = m2.decode();
    std::array<unsigned, 3> exp2{x, y, z};
    exp2[d] = (exp2[d] - k) & 0xff;
    CHECK(c2[0] == exp2[0]);
    CHECK(c2[1] == exp2[1]);
    CHECK(c2[2] == exp2[2]);
  }
}

TEST_CASE("Z-order successor/predecessor is integer ++/--") {
  using M = Morton<2, 5>;
  for (std::uint32_t i = 0; i + 1 < (1u << 10); ++i) {
    M a = M::from_code(i);
    M b = M::from_code(i + 1);
    CHECK(a < b);
    M aa = a;
    ++aa;
    CHECK(aa == b);
    M bb = b;
    --bb;
    CHECK(bb == a);
  }
  // wrap at the ends of the field
  M zero = M::from_code(0);
  M top = zero;
  --top;
  CHECK(top.code() == M::field_mask);
  ++top;
  CHECK(top == zero);
}

TEST_CASE("saturating and checked arithmetic clamp instead of wrapping") {
  using M = Morton<2, 8>;  // coord_max == 255
  SUBCASE("add_sat clamps at coord_max") {
    M m = M::encode(std::uint8_t(250), std::uint8_t(7));
    m.add_sat(0, 100);
    CHECK(m.get(0) == 255);
    CHECK(m.get(1) == 7);  // other axis untouched
  }
  SUBCASE("sub_sat clamps at 0") {
    M m = M::encode(std::uint8_t(3), std::uint8_t(7));
    m.sub_sat(0, 100);
    CHECK(m.get(0) == 0);
    CHECK(m.get(1) == 7);
  }
  SUBCASE("try_add refuses overflow and leaves code unchanged") {
    M m = M::encode(std::uint8_t(250), std::uint8_t(7));
    CHECK_FALSE(m.try_add(0, 100));
    CHECK(m.get(0) == 250);
    CHECK(m.try_add(0, 5));
    CHECK(m.get(0) == 255);
  }
  SUBCASE("try_sub refuses underflow") {
    M m = M::encode(std::uint8_t(3), std::uint8_t(7));
    CHECK_FALSE(m.try_sub(0, 100));
    CHECK(m.get(0) == 3);
    CHECK(m.try_sub(0, 3));
    CHECK(m.get(0) == 0);
  }
  SUBCASE("saturating matches wrapping when no overflow") {
    M a = M::encode(std::uint8_t(10), std::uint8_t(20));
    M b = a;
    a.add(0, 5);
    b.add_sat(0, 5);
    CHECK(a == b);
  }
}

TEST_CASE("add is consistent with repeated inc") {
  using M = Morton<2, 10>;
  std::mt19937_64 rng(3);
  for (int i = 0; i < 2000; ++i) {
    std::uint16_t x = rng() & 1023, y = rng() & 1023;
    std::uint16_t k = rng() % 50;
    M a = M::encode(x, y);
    a.add(0, k);
    M b = M::encode(x, y);
    for (int j = 0; j < k; ++j)
      b.inc(0);
    CHECK(a == b);
  }
}
