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
./build/benchmarks/morton_bench_batch             # SIMD/batch (cache vs out-of-cache)
cmake --build build --target docs                 # Doxygen -> docs/doxygen/html

# software (non-BMI2) path — must also pass:
cmake -S . -B build_nobmi2 -DMORTON_ENABLE_BMI2=OFF && cmake --build build_nobmi2 -j

# Python: proper wheel via scikit-build-core (pip install . from repo root)
pip install .                                     # builds extension, bundles .so
python -m pytest bindings/python/tests -q
# or dev loop without installing:
cmake --build build --target mortonarith_c        # drops .so into bindings/python/mortonarith/
PYTHONPATH=bindings/python python3 -m pytest bindings/python/tests -q
```

Note: the wheel build is driven by the root `pyproject.toml` → `MORTON_WHEEL=ON` in CMakeLists.txt, which builds only the binding and `return()`s early. Don't remove that early return or the wheel will also try to build tests/benchmarks.

A single doctest case can be run with `./build/tests/morton_tests --test-case="<name>"` (doctest filters, e.g. `--subcase=...`, `--list-test-cases`).

## Architecture (new core)

Header dependency order: `morton.hpp` ← {`iterate.hpp`, `batch.hpp`, `octree.hpp`}.

### `morton/morton.hpp` — `Morton<unsigned Dim, unsigned Bits>`

Interleaves `Dim` coordinates of `Bits` bits each. `Dim*Bits <= 64` uses a built-in `code_type`; **65–128 bits use `__uint128_t`** (guarded by `MORTON_HAS_INT128` / `__SIZEOF_INT128__`), so `Morton3D32` (96-bit) and `Morton2D64` (128-bit) work. Key design points:

- **Encode/decode** go through `deposit`/`extract`. For `code_bits <= 64` they use BMI2 `_pdep_u64`/`_pext_u64` under `#if defined(__BMI2__)`, else a software fallback (`detail::spread_sw`/`compact_sw`). For `code_bits > 64` the software path is always used (PDEP is 64-bit only). Paths are cross-checked in tests; the BMI2=OFF build verifies the fallback emits no pdep/pext.
- **`constexpr`**: `deposit`/`extract`/`encode`/`decode`/arithmetic are `constexpr`. They call the BMI2 intrinsic only when `!detail::is_consteval()` (which uses `__builtin_is_constant_evaluated()`, available at C++17 on GCC/Clang); in constant evaluation the software path runs. Don't "simplify" by removing the `is_consteval()` guard — the intrinsics aren't constexpr.
- **The headline arithmetic** (`add`/`sub`/`inc`/`dec`/`neighbor`, O(1), branchless) operates on one axis' interleaved bits: to add to axis `d`, fill non-axis bits with 1s (`code | ~M`) so carries ripple across the gaps, add the dilated increment, keep only axis `d` (`& M`), OR back the others. **Wraps mod `2^Bits` per axis.** `add_sat`/`sub_sat`/`try_add`/`try_sub` are the non-wrapping variants.
- **Neighbour/hierarchy helpers**: `face_neighbors()` (2·Dim, von Neumann), `all_neighbors()` (`3^Dim-1`, Moore), `ancestor(level)`/`child(level,oct)`/`child_index(level)` for octree navigation.
- `field_mask`/`coord_max` use `sizeof(T)*CHAR_BIT` (not `numeric_limits::digits`, which isn't specialised for `__int128` under strict `-std=c++17`) to avoid `1<<width` UB. Per-axis masks are an `inline static constexpr std::array axis_masks_` (constexpr-usable *and* fast at runtime).

### `morton/iterate.hpp` — region traversal + range search

- `for_each_in_box` — row-major (axis 0 fastest) odometer of O(1) `inc`/`set`; never re-encodes the full tuple. **Profiling note:** for *dense* sweeps, per-cell re-encoding is faster (independent PDEPs pipeline; the arithmetic walk is a serial dependency chain). Its value is convenience + the Z-order variant.
- `for_each_in_box_zorder` — increasing-code order via `detail::litmax_bigmin` (Tropf-Herzog). Also exposes `bigmin_in_box` / `litmax_in_box` for callers driving their own range scans. Validated exhaustively vs brute force.

### `morton/batch.hpp` — vectorised bulk ops

`batch::add/sub/step/encode2/encode3` over arrays. The masked-add loop auto-vectorises (AVX2 `vpaddq`/`vpand`/`vpor`); ~1.7× over scalar when cache-resident, memory-bound parity otherwise (`benchmarks/bench_batch.cpp`).

### `morton/octree.hpp` — `Octree<Dim,Bits,T>`

Linear octree/quadtree over `std::map<Morton,Cell>` (Cell = `{level, value}`), keyed by cell origin (low `level*Dim` bits zero). `find` uses `upper_bound` + `ancestor(level)` check; `face_neighbor` steps one cell-width with `try_add`/`try_sub` then `find`; `refine` splits into `2^Dim` children via `child()`. This is the new-core replacement for `legacy/octree.hpp`.

## Conventions / gotchas

- Everything is templated on compile-time `Dim` and `Bits`; that is what makes masks/shifts free. Don't add runtime dimension/bit-width to the hot path.
- The arithmetic advantage is for **scattered/data-dependent neighbour access**, not dense sweeps — keep this framing in docs and benchmarks. Don't market the library as a faster encoder; it is at parity with libmorton on encode/decode.
- `third_party/` vendors `doctest.h` and `libmorton/` purely for tests/benchmarks; they are not part of the shipped library.
- Python bindings are deliberately dependency-free: a C ABI shim (`bindings/morton_c.cpp`, `extern "C"`, bulk array functions) + a `ctypes` wrapper. Supported configs are `(2,32) (2,16) (3,21) (3,16)`; adding one means adding a `DEFINE_2D/3D` instantiation *and* an entry in `_CONFIG` in `__init__.py`.

## Legacy (`legacy/`)

`legacy/octree.hpp` (linear octree over a `std::map<Morton,Cell>`), `legacy/morton.hpp` (CRTP `Morton` over arbitrary-width `BitArrayBase`), `legacy/bitarray.hpp`. These require **C++20** (`<bit>` `std::countr_zero`, fold expressions) and have no build target; the octree visualization notebook is `legacy/octree_visualize.ipynb` (dumps `llx lly urx ury` rows to `temp.dat`, plots with matplotlib). Porting the octree onto the new core is a roadmap item.
