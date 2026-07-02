#include <array>
#include <cstdint>

#include "doctest.h"
#include "morton_octree/octree.hpp"

using namespace morton_octree;

TEST_CASE("uniform quadtree: point location is correct everywhere") {
  using Tree = Octree<2, 6, int>;  // 64x64 domain
  using M = Tree::key_type;
  Tree t;

  const unsigned level = 2;  // 4x4 cells
  const std::uint8_t step = 1 << level;
  int id = 0;
  for (std::uint8_t y = 0; y < 64; y += step)
    for (std::uint8_t x = 0; x < 64; x += step)
      t.insert(M::encode(x, y), level, id++);

  CHECK(t.size() == 16 * 16);

  // every point maps to the leaf whose origin is the block corner
  for (std::uint8_t px = 0; px < 64; ++px)
    for (std::uint8_t py = 0; py < 64; ++py) {
      auto it = t.find(M::encode(px, py));
      REQUIRE(it != t.end());
      auto b = t.bounds(it);
      CHECK(b[0][0] == std::uint8_t(px & ~(step - 1)));
      CHECK(b[0][1] == std::uint8_t(py & ~(step - 1)));
      CHECK(b[1][0] == std::uint8_t(b[0][0] + step - 1));
      CHECK(b[1][1] == std::uint8_t(b[0][1] + step - 1));
    }
}

TEST_CASE("face_neighbor returns the adjacent leaf") {
  using Tree = Octree<2, 6, int>;
  using M = Tree::key_type;
  Tree t;
  const unsigned level = 2;
  const std::uint8_t step = 1 << level;
  for (std::uint8_t y = 0; y < 64; y += step)
    for (std::uint8_t x = 0; x < 64; x += step)
      t.insert(M::encode(x, y), level, 0);

  auto leaf = t.find(M::encode(std::uint8_t(20), std::uint8_t(20)));
  REQUIRE(leaf != t.end());

  auto right = t.face_neighbor(leaf, 0, +1);
  REQUIRE(right != t.end());
  CHECK(t.bounds(right)[0][0] == 24);
  CHECK(t.bounds(right)[0][1] == 20);

  auto left = t.face_neighbor(leaf, 0, -1);
  REQUIRE(left != t.end());
  CHECK(t.bounds(left)[0][0] == 16);

  auto up = t.face_neighbor(leaf, 1, +1);
  REQUIRE(up != t.end());
  CHECK(t.bounds(up)[0][1] == 24);

  // edge cell has no neighbour past the domain
  auto edge = t.find(M::encode(std::uint8_t(60), std::uint8_t(60)));
  CHECK(t.face_neighbor(edge, 0, +1) == t.end());
}

TEST_CASE("refinement keeps a valid partition and point location works") {
  using Tree = Octree<2, 6, int>;
  using M = Tree::key_type;
  Tree t;

  // one root cell covering the whole 64x64 domain
  t.insert(M::encode(std::uint8_t(0), std::uint8_t(0)), 6, 0);
  CHECK(t.size() == 1);

  // refine the root into four level-5 children
  t.refine(t.begin(), [](M, unsigned oct) { return int(oct); });
  CHECK(t.size() == 4);

  // refine the child containing (40,40) once more
  auto child = t.find(M::encode(std::uint8_t(40), std::uint8_t(40)));
  REQUIRE(child != t.end());
  CHECK(child->second.level == 5);
  t.refine(child, [](M, unsigned oct) { return int(100 + oct); });
  CHECK(t.size() == 4 - 1 + 4);  // replaced one with four

  // partition still covers every point exactly once
  for (std::uint8_t px = 0; px < 64; px += 1)
    for (std::uint8_t py = 0; py < 64; py += 1)
      REQUIRE(t.find(M::encode(px, py)) != t.end());

  // (40,40) is now inside a level-4 leaf
  auto fine = t.find(M::encode(std::uint8_t(40), std::uint8_t(40)));
  CHECK(fine->second.level == 4);
}

TEST_CASE("3D octree point location") {
  using Tree = Octree<3, 5, int>;  // 32^3
  using M = Tree::key_type;
  Tree t;
  const unsigned level = 1;  // 2x2x2 cells
  const std::uint8_t step = 1 << level;
  for (std::uint8_t z = 0; z < 32; z += step)
    for (std::uint8_t y = 0; y < 32; y += step)
      for (std::uint8_t x = 0; x < 32; x += step)
        t.insert(M::encode(x, y, z), level, 0);

  for (std::uint8_t v = 0; v < 32; v += 3) {
    auto it = t.find(M::encode(v, std::uint8_t(v + 1 < 32 ? v + 1 : v), std::uint8_t(31 - v)));
    REQUIRE(it != t.end());
  }
}
