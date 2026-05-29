# morton-octree

A linear octree/quadtree built on the [morton-arithmetic](../README.md) library.

> **Status: early scaffold.** This is being split out into its own project.
> See [PLAN.md](PLAN.md) for the design and roadmap (2:1 balancing, cross-level
> neighbour queries, bulk construction, alternative storage backends).

Leaves are stored in a `std::map` keyed by their Morton origin, so the container
is already in Z-order. Point location, neighbour finding and refinement are
expressed with the morton core's O(1) arithmetic and hierarchy helpers rather
than decode/modify/encode.

```cpp
#include "morton_octree/octree.hpp"
using morton_octree::Octree;

Octree<2, 6, int> tree;                  // 64x64 quadtree of ints
tree.insert(/*coords*/{0, 0}, /*level*/6, 0);
tree.refine(tree.begin(), [](auto code, unsigned oct) { return int(oct); });
auto leaf = tree.find({40, 40});         // point location
auto right = tree.face_neighbor(leaf, /*axis*/0, /*dir*/+1);
auto box = tree.bounds(leaf);            // inclusive integer [lo, hi]
```

## Build

Depends on the installed `morton` library:

```bash
cmake -S . -B build -Dmorton_DIR=<prefix>/lib/cmake/morton
cmake --build build -j
ctest --test-dir build --output-on-failure
```

While developing inside this monorepo you can instead point at the sibling
source tree:

```bash
cmake -S . -B build -DMORTON_INCLUDE_DIR=../include
```

## License

MIT.
