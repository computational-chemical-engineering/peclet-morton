# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A header-only **C++17** library for Morton (Z-order) codes whose distinguishing feature is **arithmetic directly in Morton space** — incrementing/adding along a single axis, neighbour finding, and Z-order stepping without the usual `decode → modify → re-encode` round trip. The original prototype (an arbitrary-width `BitArray`, a wide-code `Morton`, and an octree) lives in `legacy/`; the new fast core supersedes it for codes that fit in 64 bits.

Read `docs/EVALUATION.md` for the honest assessment of what this contributes versus libmorton/morton-nd, and `docs/ROADMAP.md` for the plan.

## Build / test / benchmark

`DEVELOPMENT.md` is the canonical fresh-machine setup (prerequisites, ISA/platform notes, the full CUDA-via-clang toolkit-assembly recipe); the commands below are the day-to-day subset.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure       # or ./build/tests/morton_tests
./build/benchmarks/morton_bench                   # vs libmorton
./build/benchmarks/morton_bench_batch             # SIMD/batch (cache vs out-of-cache)
cmake --build build --target docs                 # Doxygen -> docs/doxygen/html

# software (non-BMI2) path — must also pass:
cmake -S . -B build_nobmi2 -DMORTON_ENABLE_BMI2=OFF && cmake --build build_nobmi2 -j

# runtime-dispatch path (single portable binary, no -mbmi2, uses PDEP/AVX-512
# when the CPU has them) — what the wheels build:
cmake -S . -B build_rt -DMORTON_ENABLE_BMI2=OFF -DMORTON_ENABLE_RUNTIME_DISPATCH=ON && cmake --build build_rt -j

# AVX-512 has no hardware on most dev machines / CI; validate it under Intel SDE
# (userspace emulator, runs on AMD too). -skx => AVX-512 batch path executes:
sde64 -skx -- ./build/tests/morton_tests
# runtime-dispatch build across emulated CPUs (BMI2 present/absent):
sde64 -snb -- ./build_rt/tests/morton_tests   # no BMI2/AVX-512 -> software fallback
sde64 -hsw -- ./build_rt/tests/morton_tests   # BMI2, no AVX-512 -> PDEP path
sde64 -skx -- ./build_rt/tests/morton_tests   # BMI2 + AVX-512

# portable Kokkos GPU backend (opt-in; needs a Kokkos prefix on CMAKE_PREFIX_PATH
# from ../tools/bootstrap_deps.sh <nvidia-cuda|host-openmp|lumi-hip>):
cmake -S . -B build_kokkos -DMORTON_ENABLE_KOKKOS=ON \
      -DCMAKE_PREFIX_PATH="$PWD/../extern/install/host-openmp" && cmake --build build_kokkos -j
ctest --test-dir build_kokkos -R morton_kokkos_tests --output-on-failure

# Python: proper wheel via scikit-build-core (pip install . from repo root)
pip install .                                     # builds extension, bundles .so
python -m pytest bindings/python/tests -q
# or dev loop without installing:
cmake --build build --target mortonarith_c        # drops .so into bindings/python/peclet/morton/
PYTHONPATH=bindings/python python3 -m pytest bindings/python/tests -q
```

Note: the wheel build is driven by the root `pyproject.toml` → `MORTON_WHEEL=ON` in CMakeLists.txt, which builds only the binding and `return()`s early. Don't remove that early return or the wheel will also try to build tests/benchmarks.

A single doctest case can be run with `./build/tests/morton_tests --test-case="<name>"` (doctest filters, e.g. `--subcase=...`, `--list-test-cases`).

## Architecture (new core)

Header dependency order: `morton.hpp` ← {`iterate.hpp`, `simd.hpp` ← `batch.hpp`}. (The octree was split
out into the sibling `octree/` subproject — see below — and is no longer a core header.)

### `morton/morton.hpp` — `Morton<unsigned Dim, unsigned Bits>`

Interleaves `Dim` coordinates of `Bits` bits each. `Dim*Bits <= 64` uses a built-in `code_type`; **65–128 bits use `__uint128_t`** (guarded by `MORTON_HAS_INT128` / `__SIZEOF_INT128__`), so `Morton3D32` (96-bit) and `Morton2D64` (128-bit) work. Key design points:

- **Encode/decode** go through `deposit`/`extract`. For `code_bits <= 64` they use BMI2 `_pdep_u64`/`_pext_u64` under `#if defined(__BMI2__)`, else a software fallback (`detail::spread_sw`/`compact_sw`). For `code_bits > 64` the software path is always used (PDEP is 64-bit only). Paths are cross-checked in tests; the BMI2=OFF build verifies the fallback emits no pdep/pext.
- **Runtime dispatch** (`MORTON_ENABLE_RUNTIME_DISPATCH`, opt-in via the CMake option / compile define; gated to `MORTON_X86_RUNTIME_DISPATCH` on x86-64 GCC/Clang when `-mbmi2` is *not* set): a third `deposit`/`extract` branch calls `detail::pdep_u64_hw`/`pext_u64_hw` — `__attribute__((target("bmi2")))` helpers, so PDEP is emitted without a global `-mbmi2` — guarded by the cached `detail::cpu_has_bmi2()` (`__builtin_cpu_supports`). This yields one portable binary that's BMI2-fast at runtime. **Don't enable it by default for the plain BMI2=OFF build** — that build is contractually pdep-free (a test greps for it). It self-disables under `-mbmi2`/CUDA/non-x86. `MORTON_X86` / `cpu_has_bmi2`/`cpu_has_avx2`/`cpu_has_avx512f` live in `detail` and are shared with the SIMD path.
- **`constexpr`**: `deposit`/`extract`/`encode`/`decode`/arithmetic are `constexpr`. They call the BMI2 intrinsic only when `!detail::is_consteval()` (which uses `__builtin_is_constant_evaluated()`, available at C++17 on GCC/Clang); in constant evaluation the software path runs. Don't "simplify" by removing the `is_consteval()` guard — the intrinsics aren't constexpr.
- **The headline arithmetic** (`add`/`sub`/`inc`/`dec`/`neighbor`, O(1), branchless) operates on one axis' interleaved bits: to add to axis `d`, fill non-axis bits with 1s (`code | ~M`) so carries ripple across the gaps, add the dilated increment, keep only axis `d` (`& M`), OR back the others. **Wraps mod `2^Bits` per axis.** `add_sat`/`sub_sat`/`try_add`/`try_sub` are the non-wrapping variants.
- **Neighbour/hierarchy helpers**: `face_neighbors()` (2·Dim, von Neumann), `all_neighbors()` (`3^Dim-1`, Moore), `ancestor(level)`/`child(level,oct)`/`child_index(level)` for octree navigation.
- `field_mask`/`coord_max` use `sizeof(T)*CHAR_BIT` (not `numeric_limits::digits`, which isn't specialised for `__int128` under strict `-std=c++17`) to avoid `1<<width` UB. Per-axis masks are an `inline static constexpr std::array axis_masks_` (constexpr-usable *and* fast at runtime).

### `morton/iterate.hpp` — region traversal + range search

- `for_each_in_box` — row-major (axis 0 fastest) odometer of O(1) `inc`/`set`; never re-encodes the full tuple. **Profiling note:** for *dense* sweeps, per-cell re-encoding is faster (independent PDEPs pipeline; the arithmetic walk is a serial dependency chain). Its value is convenience + the Z-order variant.
- `for_each_in_box_zorder` — increasing-code order via `detail::litmax_bigmin` (Tropf-Herzog). Also exposes `bigmin_in_box` / `litmax_in_box` for callers driving their own range scans. Validated exhaustively vs brute force.

### `morton/batch.hpp` — vectorised bulk ops

`batch::add/sub/step/encode2/encode3/decode2/decode3` over arrays. The masked-add loop auto-vectorises (AVX2 `vpaddq`/`vpand`/`vpor`); ~1.7× over scalar when cache-resident, memory-bound parity otherwise (`benchmarks/bench_batch.cpp`). On x86-64 GCC/Clang (`MORTON_X86`) each function **runtime-dispatches** to the explicit AVX-512 kernels in `simd.hpp` when `detail::cpu_has_avx512f()` and the code is 64-bit — transparently, with bit-for-bit identical results — else it falls through to this scalar loop.

### `morton/simd.hpp` — explicit AVX-512 batch kernels

Hand-written AVX-512 "magic-bits" spread/compact (`detail::avx512::{encode2,decode2,encode3,decode3,add64,sub64}`), 8 codes per iteration, with a scalar tail (so `n` need not be a multiple of 8). Each kernel is `__attribute__((target("avx512f")))`, so the TU needs **no** global `-mavx512f` and the rest of the binary stays portable; they only run when `cpu_has_avx512f()` is true. Coverage is the 64-bit-code, 32-bit-coord layouts: `encode/decode` for `(2,32)` and `(3,21)`, and `add/sub` for any 64-bit code (`(2,32)`,`(3,21)`,`(3,16)`); other layouts use the scalar path. **Validate changes here under Intel SDE** (`sde64 -skx`) — there's usually no AVX-512 hardware locally, and the tests compare the SIMD output against the scalar reference exhaustively.

### `morton/wide_uint.hpp` — codes wider than 128 bits

`detail::wide_uint<W>` is a minimal fixed-width (`W` × `u64`, little-endian) unsigned providing exactly the operators `Morton` uses (`+ - & | ^ ~ << >>`, comparisons, conversions to `u64`/`__int128`). `uint_for` selects it automatically when `Dim*Bits` exceeds the builtin width (64, or 128 with `__int128`). Cap is `MORTON_MAX_BITS` (default 256). The Morton class body is unchanged — it's all operators — so 192-bit (`Morton<3,64>`) and 256-bit (`Morton<2,128>`) "just work". Don't tune `wide_uint` for speed unless it becomes a hot path; it's deliberately simple.

### `include/morton/kokkos.hpp` — portable Kokkos GPU backend (shares the core via `MORTON_HD`)

The core's functions are prefixed with `MORTON_HD`. It now resolves to **`KOKKOS_FUNCTION`** when `MORTON_ENABLE_KOKKOS` is set and Kokkos is in the translation unit (so device/host attributes track whichever Kokkos backend — CUDA/HIP/OpenMP/Serial — is active), and falls back to `__host__ __device__` under `__CUDACC__`/`__HIPCC__`, and to nothing for an ordinary host build (so the CPU library is byte-for-byte unchanged). This lets the Kokkos kernels call the **same** `Morton<Dim,Bits>` code — no second implementation. Two gotchas baked into the core (still required under the Kokkos CUDA/HIP device pass):
- `deposit`/`extract` guard PDEP/PEXT with `#if defined(__BMI2__) && !defined(__CUDA_ARCH__)` — never emit the x86 intrinsic in device code, even when the host is compiled `-mbmi2`. Don't drop the `!defined(__CUDA_ARCH__)`.
- `axis_mask(d)` computes the mask under `__CUDA_ARCH__` instead of reading the host `static constexpr axis_masks_` array (can't reference a host static from device).

`morton::kokkos` offers `Kokkos::View`-based bulk ops — `encode2/3`, `decode2/3`, per-axis `add`/`sub`/`step` (execution space deduced from the View) — plus `*_host` raw-pointer convenience wrappers (stage to device Views, run, copy back). `tests/kokkos/test_kokkos.cpp` validates device output bit-for-bit against the scalar `Morton<>` reference for all four binding layouts `(2,32) (2,16) (3,21) (3,16)`; `benchmarks/bench_kokkos.cpp` shows device-resident ~51 GMops/s (CUDA) vs PCIe-bound host round-trip.

Build: opt-in with `-DMORTON_ENABLE_KOKKOS=ON` + `find_package(Kokkos CONFIG)` against the suite's bootstrapped prefix (`extern/install/<backend>` from `tools/bootstrap_deps.sh`) on `CMAKE_PREFIX_PATH`, exactly like `flow`/`dem`. Device sources stay plain `.cpp` — Kokkos 5.x routes them through the launch compiler (see `../cmake/SuiteKokkos.cmake`). The plain (non-Kokkos) build never touches Kokkos and is unchanged. **The raw-CUDA backend was retired** (was `cuda/`, `morton::cuda`); the last raw-CUDA tree is at the `pre-cuda-retirement` git tag.

### Octree → sibling `octree/` project

The octree is **no longer part of this library**. It moved to `octree/` (`morton_octree::Octree`, `octree/include/morton_octree/octree.hpp`), a separate project that depends on `morton::morton`. See `octree/PLAN.md`. `legacy/octree.hpp` remains the original prototype/reference.

## Packaging / distribution

- **CMake**: `install` exports a package config (`cmake/morton-config.cmake.in` → `find_package(morton CONFIG)` → `morton::morton`). The `-mbmi2` interface flag is guarded by `$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,AppleClang>` so it's safe to propagate to consumers.
- **Conan** `conanfile.py`, **vcpkg** `packaging/vcpkg/morton/` (set REF/SHA512 at release).
- **Python wheels**: `.github/workflows/release.yml` + cibuildwheel build with `MORTON_ENABLE_BMI2=OFF -DMORTON_ENABLE_RUNTIME_DISPATCH=ON` (via `CMAKE_ARGS`) — no global `-mbmi2` (so no SIGILL on old CPUs) but PDEP/PEXT and the AVX-512 batch kernels still kick in at runtime via `target`-attribute helpers + CPUID. One wheel, portable *and* fast. (On MSVC/non-x86 the dispatch self-disables → software path.) Source `pip install .` stays BMI2-on. The C ABI shim (`bindings/morton_c.cpp`) routes 64-bit-code configs through `batch::` so Python gets the dispatched paths; `(2,16)` (32-bit code) keeps the per-element path.

## Conventions / gotchas

- Everything is templated on compile-time `Dim` and `Bits`; that is what makes masks/shifts free. Don't add runtime dimension/bit-width to the hot path.
- The arithmetic advantage is for **scattered/data-dependent neighbour access**, not dense sweeps — keep this framing in docs and benchmarks. Don't market the library as a faster encoder; it is at parity with libmorton on encode/decode.
- `third_party/` vendors `doctest.h` and `libmorton/` purely for tests/benchmarks; they are not part of the shipped library.
- Python bindings are deliberately dependency-free: a C ABI shim (`bindings/morton_c.cpp`, `extern "C"`, bulk array functions) + a `ctypes` wrapper. Supported configs are `(2,32) (2,16) (3,21) (3,16)`; adding one means adding a `DEFINE_2D/3D` instantiation *and* an entry in `_CONFIG` in `__init__.py`.

## Legacy (`legacy/`)

`legacy/octree.hpp` (linear octree over a `std::map<Morton,Cell>`), `legacy/morton.hpp` (CRTP `Morton` over arbitrary-width `BitArrayBase`), `legacy/bitarray.hpp`. These require **C++20** (`<bit>` `std::countr_zero`, fold expressions) and have no build target; the octree visualization notebook is `legacy/octree_visualize.ipynb` (dumps `llx lly urx ury` rows to `temp.dat`, plots with matplotlib). Porting the octree onto the new core is a roadmap item.
