#include "doctest.h"
#include "morton/morton.hpp"

using namespace morton;

// Everything below is evaluated at compile time: if the software constexpr
// path were broken these static_asserts would fail to compile.

constexpr auto c_enc = Morton<3, 21>::encode(100u, 200u, 300u);
static_assert(c_enc.get(0) == 100, "");
static_assert(c_enc.get(1) == 200, "");
static_assert(c_enc.get(2) == 300, "");

constexpr auto c_add = [] {
  auto m = Morton<2, 16>::encode(std::uint16_t(1000), std::uint16_t(2000));
  m.add(0, 234);
  m.sub(1, 500);
  return m;
}();
static_assert(c_add.get(0) == 1234, "");
static_assert(c_add.get(1) == 1500, "");

constexpr auto c_round = Morton<2, 8>::encode(std::uint8_t(0xAB), std::uint8_t(0xCD));
static_assert(c_round.decode()[0] == 0xAB, "");
static_assert(c_round.decode()[1] == 0xCD, "");

// A compile-time lookup table built from the library (the motivating use case).
constexpr std::array<std::uint32_t, 8> row0 = [] {
  std::array<std::uint32_t, 8> t{};
  auto m = Morton<2, 8>::encode(std::uint8_t(0), std::uint8_t(0));
  for (int x = 0; x < 8; ++x) {
    t[x] = std::uint32_t(m.code());
    m.inc(0);
  }
  return t;
}();
static_assert(row0[0] == 0 && row0[1] == 1 && row0[2] == 4 && row0[3] == 5, "");

TEST_CASE("constexpr results agree with runtime results") {
  // Runtime evaluation (BMI2 path) must match the compile-time values above.
  auto m = Morton<3, 21>::encode(100u, 200u, 300u);
  CHECK(m.code() == c_enc.code());

  auto t = Morton<2, 16>::encode(std::uint16_t(1000), std::uint16_t(2000));
  t.add(0, 234);
  t.sub(1, 500);
  CHECK(t.code() == c_add.code());
}
