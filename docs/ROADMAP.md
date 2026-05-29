# Roadmap / plan

Where this library stands and what to do with it.

## Status today (v0.2)

The v0.1 core plus the entire near/medium-term roadmap below is now
**implemented and tested** (27 doctest cases, ~1.4M assertions; pytest; CI; docs):

- Header-only C++17 core (`Morton<Dim, Bits>`): BMI2 encode/decode with portable
  software fallback; O(1) per-axis `inc/dec/add/sub`, `neighbor`, `set`, Z-order
  `++`/`--`, comparisons.
- ✅ **(1) CI matrix** — `.github/workflows/ci.yml`: gcc/clang/MSVC ×
  {BMI2 on/off} × {Debug/Release} running `ctest`, plus a Python-version matrix
  and a Doxygen build. The software fallback is exercised on the BMI2=OFF cells.
- ✅ **(2) `constexpr` encode/decode/arithmetic** — selects the software path at
  compile time via `__builtin_is_constant_evaluated()` (works at C++17); see
  `tests/test_constexpr.cpp` for compile-time lookup-table construction.
- ✅ **(3) Saturating / checked arithmetic** — `add_sat`, `sub_sat`, `try_add`,
  `try_sub` clamp or refuse instead of wrapping.
- ✅ **(4) Real Python build** — `scikit-build-core`; `pip install .` builds the
  extension and produces a proper wheel (`mortonarith-*.whl`) with the `.so`
  bundled. (PyPI publishing still pending — see below.)
- ✅ **(5) Docs site** — Doxygen config (`docs/Doxyfile`) + `docs` CMake target.
- ✅ **(6) Wider codes (> 64 bits)** — `__uint128_t` storage where available, so
  3D 32-bit (96-bit code) and 2D 64-bit (128-bit code) work; software
  encode/decode for >64 bits, native 128-bit arithmetic. Aliases `Morton3D32`,
  `Morton2D64`.
- ✅ **(7) LITMAX companion to BIGMIN** — `morton/iterate.hpp` exposes the full
  Tropf-Herzog pair (`bigmin_in_box`, `litmax_in_box`, `litmax_bigmin`),
  validated exhaustively against brute force.
- ✅ **(8) Neighbour-set + hierarchy helpers** — `face_neighbors()` (von
  Neumann), `all_neighbors()` (Moore, `3^Dim-1`), `ancestor`/`child`/
  `child_index`.
- ✅ **(9) Octree on the new core** — `morton/octree.hpp`: a linear octree/
  quadtree (`std::map` keyed by Morton origin) with point location, face
  neighbours and refinement expressed via the core arithmetic.
- ✅ **(10) SIMD batch arithmetic** — `morton/batch.hpp` (`add`/`sub`/`step`/
  `encode2`/`encode3`) auto-vectorises (AVX2 `vpaddq`/`vpand`/`vpor`). Profiling
  (`benchmarks/bench_batch.cpp`) confirms ~1.7× over scalar when cache-resident
  and memory-bound parity out of cache — exactly the predicted behaviour.

## Remaining / future work

- **Publish to PyPI** with `cibuildwheel` to produce manylinux/macos/windows
  wheels (the build is wheel-ready; only the release pipeline is missing).
- **C++ packaging**: vcpkg / Conan / CMake package config export so
  `find_package(morton)` works from an install tree (install rules exist).
- **Explicit AVX-512** batch encode (libmorton-style) for the cases where the
  bulk path is *not* memory-bound (cache-resident transforms).
- **Codes > 128 bits** via a small fixed word-array backend (the arithmetic
  generalises with inter-word carry); reuse the legacy `BitArray` only as the
  slow reference.
- **Hilbert curve** option alongside Morton (better locality); the arithmetic
  framework (BIGMIN/LITMAX, neighbour stepping) largely carries over.
- **GPU / SYCL** batch kernels for the encode/arithmetic loops.
- **2:1 balancing** and parent/child iteration on the octree (the legacy
  prototype's unfinished `balanceTree`).

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
