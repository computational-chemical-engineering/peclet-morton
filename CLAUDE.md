# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A header-only **C++17** library for Morton (Z-order) codes whose distinguishing feature is **arithmetic directly in Morton space** — incrementing/adding along a single axis, neighbour finding, and Z-order stepping without the usual `decode → modify → re-encode` round trip. The original prototype (an arbitrary-width `BitArray`, a wide-code `Morton`, and an octree) lives in `legacy/`; the new fast core supersedes it for codes that fit in 64 bits.

Read `docs/EVALUATION.md` for the honest assessment of what this contributes versus libmorton/morton-nd, and `docs/ROADMAP.md` for the plan.

## Build / test / benchmark

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure       # or ./build/tests/morton_tests
./build/benchmarks/morton_bench                   # vs libmorton

# software (non-BMI2) path — must also pass:
cmake -S . -B build_nobmi2 -DMORTON_ENABLE_BMI2=OFF && cmake --build build_nobmi2 -j

# Python bindings (ctypes + NumPy, no pybind11):
cmake --build build --target mortonarith_c        # drops .so into bindings/python/mortonarith/
PYTHONPATH=bindings/python python3 -m pytest bindings/python/tests -q
```

A single doctest case can be run with `./build/tests/morton_tests --test-case="<name>"` (doctest filters, e.g. `--subcase=...`, `--list-test-cases`).

## Architecture (new core)

Dependency order: `include/morton/morton.hpp` ← `include/morton/iterate.hpp`.

### `morton/morton.hpp` — `Morton<unsigned Dim, unsigned Bits>`

Interleaves `Dim` coordinates of `Bits` bits each; `static_assert(Dim*Bits <= 64)`. The code is stored in the smallest unsigned int that fits (`code_type`), coordinates in `coord_type`. Key design points:

- **Encode/decode** go through `deposit`/`extract`, which use BMI2 `_pdep_u64`/`_pext_u64` under `#if defined(__BMI2__)` and a portable constexpr software fallback (`detail::spread_sw`/`compact_sw`) otherwise. The two paths are cross-checked in the tests; the `build_nobmi2` build verifies the fallback emits no pdep/pext.
- **The headline arithmetic** (`add`/`sub`/`inc`/`dec`/`neighbor`, all O(1), branchless) operates on one axis' interleaved bits. The trick: to add to axis `d`, fill the *non-axis* bits with 1s (`code | ~M`) so carries ripple across the gaps, add the dilated increment, then keep only axis `d`'s result (`& M`) and OR back the other axes. Subtraction uses `code & M` (non-axis bits 0) so borrows ripple. **Wrapping is modulo `2^Bits` per axis** by design.
- `field_mask`/`coord_max` are computed to avoid `1 << 64` UB when `Dim*Bits == 64`.
- Z-order successor/predecessor are just `++`/`--` on the integer code (the code *is* the Z-order index).
- Per-axis masks are memoised in `mask_table()` (function-local static).

### `morton/iterate.hpp` — region traversal

- `for_each_in_box` — row-major (axis 0 fastest) using an odometer of O(1) `inc`/`set` calls; never re-encodes the full tuple. **Note from profiling:** for *dense* sweeps, per-cell re-encoding is actually faster (independent PDEPs pipeline; the arithmetic walk is a serial dependency chain). This iterator's real value is convenience + the Z-order variant, not dense-fill speed.
- `for_each_in_box_zorder` — visits the box in increasing code order using `detail::bigmin` (Tropf-Herzog BIGMIN). Validated against brute force for 2D and 3D.

## Conventions / gotchas

- Everything is templated on compile-time `Dim` and `Bits`; that is what makes masks/shifts free. Don't add runtime dimension/bit-width to the hot path.
- The arithmetic advantage is for **scattered/data-dependent neighbour access**, not dense sweeps — keep this framing in docs and benchmarks. Don't market the library as a faster encoder; it is at parity with libmorton on encode/decode.
- `third_party/` vendors `doctest.h` and `libmorton/` purely for tests/benchmarks; they are not part of the shipped library.
- Python bindings are deliberately dependency-free: a C ABI shim (`bindings/morton_c.cpp`, `extern "C"`, bulk array functions) + a `ctypes` wrapper. Supported configs are `(2,32) (2,16) (3,21) (3,16)`; adding one means adding a `DEFINE_2D/3D` instantiation *and* an entry in `_CONFIG` in `__init__.py`.

## Legacy (`legacy/`)

`legacy/octree.hpp` (linear octree over a `std::map<Morton,Cell>`), `legacy/morton.hpp` (CRTP `Morton` over arbitrary-width `BitArrayBase`), `legacy/bitarray.hpp`. These require **C++20** (`<bit>` `std::countr_zero`, fold expressions) and have no build target; the octree visualization notebook is `legacy/octree_visualize.ipynb` (dumps `llx lly urx ury` rows to `temp.dat`, plots with matplotlib). Porting the octree onto the new core is a roadmap item.
