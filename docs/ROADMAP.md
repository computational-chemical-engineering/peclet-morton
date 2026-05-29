# Roadmap / plan

Where this library stands and what to do with it.

## Status today (v0.3)

The v0.1 core plus the entire near/medium-term roadmap below is now
**implemented and tested** (26 doctest cases, ~1.4M assertions; pytest; CI; docs);
the v0.2→v0.3 additions (wide codes, packaging, release pipeline, octree split)
are listed in "Done since v0.2" further down:

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

## Done since v0.2

- ✅ **Codes > 128 bits** — `morton/wide_uint.hpp` provides a fixed-width
  word-array unsigned (`+ - & | ^ ~ << >>`, comparisons) so `Morton<Dim,Bits>`
  works unchanged up to `MORTON_MAX_BITS` (default 256; raise to go wider).
  Tested for 192-bit (`Morton<3,64>`) and 256-bit (`Morton<2,128>`) codes on
  both the BMI2 and software paths.
- ✅ **C++ packaging** — CMake package-config export (`find_package(morton)` →
  `morton::morton`), verified by an external consumer build. Conan recipe
  (`conanfile.py`) and a vcpkg port (`packaging/vcpkg/morton/`).
- ✅ **PyPI release pipeline** — `cibuildwheel` config + `.github/workflows/
  release.yml` (manylinux/musllinux/macOS/Windows, trusted publishing on tag).
  **Released wheels build the portable software path** (`MORTON_ENABLE_BMI2=OFF`)
  so they never SIGILL on a non-BMI2 CPU; source installs stay BMI2-fast.
- ✅ **Octree split out** — moved to the sibling `octree/` project
  (`morton_octree::Octree`), which depends on this library. 2:1 balancing,
  cross-level neighbour queries and bulk construction are tracked in
  [`../octree/PLAN.md`](../octree/PLAN.md), not here.
- ✅ **GPU (CUDA) backend** — `cuda/` (`morton::cuda`). The core is marked
  `__host__ __device__` (`MORTON_HD`), so kernels reuse the exact CPU code path
  (PDEP guarded out of device code). encode/decode (2D/3D) + per-axis arithmetic;
  validated bit-for-bit against the CPU library on the GPU. ~51 GMops/s
  device-resident 2D-32 encode on an RTX 5080 (~33× one CPU core); one-shot host
  calls are PCIe-bound (documented).

## Done since v0.3

- ✅ **Runtime BMI2 dispatch** — `MORTON_ENABLE_RUNTIME_DISPATCH` (CMake option /
  compile define) builds a single binary *without* `-mbmi2` that still uses
  PDEP/PEXT when the running CPU has BMI2 (software fallback otherwise), via
  per-function `__attribute__((target("bmi2")))` helpers + a cached
  `__builtin_cpu_supports` check (`morton/morton.hpp`). The redistributable
  wheels now build with it on (`pyproject.toml`), so one portable wheel is both
  SIGILL-safe on old CPUs and BMI2/AVX-512-fast on new ones — removing the old
  portability/speed trade-off. Self-disables when `-mbmi2` is already set, on
  non-x86, and under CUDA.
- ✅ **Explicit AVX-512 batch encode/decode** — `morton/simd.hpp`: hand-written
  AVX-512 "magic-bits" kernels (8 codes/iteration) for `encode2`/`decode2`
  (2,32), `encode3`/`decode3` (3,21), and per-axis `add`/`sub` on any 64-bit
  code. `morton/batch.hpp` dispatches to them at runtime when the CPU has
  AVX-512F (else the auto-vectorised scalar path); they use a `target` attribute
  so no global `-mavx512f` is needed. The C ABI / NumPy bindings route through
  these, so Python gets the speedup too. **No longer ships untested**: CI runs
  the suite under Intel SDE (`-skx`) so the AVX-512 path is checked bit-for-bit
  against the scalar reference on ordinary (non-AVX-512) runners; the
  runtime-dispatch build is also exercised under `-snb`/`-hsw`/`-skx`.

## Remaining / future work

- **AVX-512 on real silicon**: SDE validates *correctness*; a benchmark on an
  actual AVX-512 host (or a GPU/cloud runner) is still wanted to confirm the
  expected throughput win for cache-resident bulk encode (`bench_batch` already
  prints the active path and a `batch::encode2` vs scalar comparison).
- **Wider SIMD configs**: the AVX-512 encode kernels currently cover the
  64-bit-code, 32-bit-coord layouts (2,32)/(3,21); 16-bit-coord configs and
  AVX2-only explicit kernels could be added if profiling warrants.
- **MSVC runtime dispatch**: the `target`-attribute + `__builtin_cpu_supports`
  mechanism is GCC/Clang; an `__cpuid`-based path would extend it to MSVC.
- **GPU follow-ups**: a Z-order **radix sort** (the usual reason to be on the
  GPU), a NumPy-device Python entry point, SYCL/HIP portability, and a GPU CI
  runner. (The CUDA encode/decode/arithmetic backend itself is now done — see
  above.) Design notes in [`HILBERT_GPU_NOTES.md`](HILBERT_GPU_NOTES.md).
- **Hilbert curve** option alongside Morton: a larger, separable effort — design
  notes in [`HILBERT_GPU_NOTES.md`](HILBERT_GPU_NOTES.md).
- **Codes > 256 bits**: just raise `MORTON_MAX_BITS`; revisit `wide_uint`
  performance (it is intentionally simple, not tuned) if that becomes a hot path.

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
