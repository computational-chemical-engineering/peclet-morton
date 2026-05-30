# morton-octree — project plan

This directory is the seed of a **separate project**, split out from the
morton-arithmetic library. The morton library should stay a small, focused
primitive (encode/decode + Z-order arithmetic + range search); a full octree is
a different concern with its own data-structure and numerics decisions, so it
lives — and will eventually ship — on its own.

It **depends on** the standalone morton library (`morton::Morton`) and uses its
O(1) arithmetic (`inc`/`add`/`neighbor`), hierarchy helpers
(`ancestor`/`child`/`child_index`) and range queries (BIGMIN/LITMAX).

## What exists today (scaffold)

`include/morton_octree/octree.hpp` — a linear octree/quadtree:
`std::map<Morton, Cell>` keyed by cell origin, with point location
(`find`), `face_neighbor`, `bounds`, and `refine`. Tested in `tests/`
(4 cases). This is enough to validate that the morton primitives are the right
substrate; it is **not** yet a finished octree.

## Roadmap for the standalone project

1. **Own identity.** Repository, namespace `morton_octree`, CMake package
   (`morton_octree::morton_octree`) with a `find_package(morton)` dependency,
   semantic versioning. (CMakeLists.txt here is the starting point.)
2. **2:1 balancing.** The original prototype's unfinished `balanceTree`: enforce
   that face-adjacent leaves differ by at most one level, refining as needed.
   This is the main numerical feature an AMR/mesh user needs. Neighbour-of-equal
   -size lookups are already O(1) via the morton arithmetic.
3. **Full neighbour queries across levels.** `face_neighbors` that return *all*
   same-or-finer neighbours across a face (not just the containing leaf), built
   on `child`/`ancestor` and the range-query helpers.
4. **Payload ergonomics.** Iterators over leaves in Z-order (free — the map is
   already Z-ordered), parent/child traversal, and a visitor API. Optional
   value-less "set" specialisation.
5. **Bulk construction.** Build from a point cloud or a predicate by sorting
   Morton codes and running a linear bottom-up merge (classic linear-octree
   construction), rather than repeated `insert`.
6. **Storage alternatives.** Benchmark `std::map` vs a sorted `std::vector` of
   `(code, cell)` (better cache behaviour, `std::lower_bound` lookups) vs a hash
   map. The arithmetic layer is storage-agnostic.
7. **Serialization / interop.** Dump leaves as `(origin, level, value)` or
   bounding boxes (the legacy `temp.dat` + matplotlib notebook workflow lives in
   `../legacy/octree_visualize.ipynb`).
8. **Tests & docs** mirroring the morton library's standards (brute-force
   validation, CI, Doxygen).

## Why separate

- Keeps the morton library dependency-free and easy to audit/adopt.
- Octree design choices (balancing policy, storage backend, payload model) are
  opinionated and shouldn't be forced on users who only want Z-order codes.
- Lets the octree iterate on its own release cadence while pinning a morton
  version.

## Relationship to `../legacy/`

`../legacy/octree.hpp` is the original arbitrary-width-BitArray octree prototype
(with the half-finished `balanceTree` and face detection). It is the reference
for features to port; this project re-implements them on the fast fixed-width
morton core.
