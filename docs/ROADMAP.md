# Roadmap / plan

Where this library stands and what to do with it. Ordered by value-to-effort.

## Status today (v0.1)

Done in this iteration:

- Header-only C++17 core (`Morton<Dim, Bits>`, `Dim*Bits ≤ 64`): BMI2
  encode/decode with portable software fallback; O(1) per-axis
  `inc/dec/add/sub`, `neighbor`, `set`, Z-order `++`/`--`, comparisons.
- Region iteration (`for_each_in_box`, `for_each_in_box_zorder` via BIGMIN).
- doctest suite (>1M assertions): cross-checked against a bit-by-bit reference
  and against libmorton; arithmetic and iterators validated against brute force.
- CMake (interface target `morton::morton`, `ctest`, install rules).
- C ABI shim + vectorised NumPy/ctypes Python package + pytest.
- Benchmarks vs libmorton; honest profiling write-up in `EVALUATION.md`.

## Near term — make it a credible open-source library

1. **CI matrix.** GitHub Actions: gcc/clang/MSVC × {BMI2 on, BMI2 off} ×
   {Debug, Release}, run `ctest` + pytest. The software fallback path must be
   exercised on a non-BMI2 build (the tests already pass there locally).
2. **`constexpr` encode/decode.** The software path is already constexpr-friendly;
   make `encode`/`decode`/arithmetic usable in constant expressions (the BMI2
   intrinsics are not constexpr, so select the software path inside
   `if (std::is_constant_evaluated())`, C++20). Lets users build compile-time
   lookup tables.
3. **Signed / saturating arithmetic options.** Currently axis arithmetic wraps
   mod 2^Bits. Add `add_saturating` and bounds-checked variants for grid code
   that must not wrap.
4. **A real Python build.** Replace the manual CMake-drops-an-`.so` step with
   `scikit-build-core` so `pip install .` produces a proper wheel; publish to
   PyPI. Consider `pybind11`/`nanobind` if a richer object API is wanted, but the
   ctypes+bulk-array design is deliberately dependency-free and fast — keep it
   unless there's a reason not to.
5. **Docs site.** The README + EVALUATION are enough to start; add Doxygen or a
   short mdBook with the neighbour/stencil examples.

## Medium term — features that justify the "arithmetic" name

6. **Wider codes (> 64 bits) without losing speed.** Reuse the legacy
   `BitArray` only as a fallback; provide a `Morton<Dim, Bits>` specialisation
   backed by `__uint128_t` (and an array of words beyond that) so 3D 32-bit and
   2D 64-bit work. The arithmetic generalises directly (carry across words).
7. **LITMAX companion to BIGMIN.** Expose the full Tropf-Herzog pair so callers
   can drive their own range-search / database-style scans, not just the
   bundled iterator.
8. **Neighbour-set helpers.** `face_neighbors()`, `all_neighbors()` (Moore /
   von Neumann), and parent/child navigation — the operations an octree actually
   calls. This is the bridge to item 9.
9. **Port the legacy octree onto the new core.** `legacy/octree.hpp` is the
   original motivating application. Re-implement it over `Morton<Dim,Bits>`; the
   `findCell`/neighbour logic becomes the fast arithmetic ops here. This turns
   the library from "primitives" into "primitives + a flagship user."
10. **SIMD batch arithmetic.** libmorton has AVX-512 encode; the per-axis
    arithmetic vectorises cleanly too (it is just masked adds). A batched
    `shift_many` could widen the Python advantage further. Profile first — the
    scalar bulk path is already memory-bound in Python.

## Positioning

Market it precisely (see `EVALUATION.md §5`): **not** a faster encoder, but the
Z-order *navigation* layer that libmorton/morton-nd intentionally omit. Realistic
homes for the work:

- Publish as a standalone small library (header-only C++ + PyPI wheel) aimed at
  octree/quadtree, stencil-on-SFC, and spatial-index authors.
- Or contribute the arithmetic + range-query operations *to* an existing library
  (e.g. as an extension layer over libmorton) so users get encode/decode and
  navigation from one place. This avoids fragmenting the ecosystem and is
  probably the higher-impact path.

## Non-goals

- Beating libmorton at raw encode/decode (it is already optimal on this
  hardware; we match it).
- Arbitrary runtime dimensions/bit-widths in the hot path — compile-time
  `Dim`/`Bits` is what makes the masks and shifts free.
