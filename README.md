# morton-arithmetic

[![PyPI version](https://img.shields.io/pypi/v/peclet-morton.svg)](https://pypi.org/project/peclet-morton/)
[![Python versions](https://img.shields.io/pypi/pyversions/peclet-morton.svg)](https://pypi.org/project/peclet-morton/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![CI](https://github.com/computational-chemical-engineering/peclet-morton/actions/workflows/ci.yml/badge.svg)](https://github.com/computational-chemical-engineering/peclet-morton/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21132433.svg)](https://doi.org/10.5281/zenodo.21132433)

A small, header-only **C++17** library for Morton (Z-order) codes whose
distinguishing feature is **arithmetic directly in Morton space**.

Most Morton libraries only *encode* (interleave coordinates into a code) and
*decode* (the inverse). This one adds the operations you actually want when you
*navigate* a Z-ordered grid:

- increment / decrement / add along a single axis,
- find a neighbour,
- step to the next or previous cell in Z-order,
- sweep an axis-aligned region,

…all **without** the usual `decode → modify → re-encode` round trip. Each is a
handful of branchless integer instructions.

```cpp
#include "morton/morton.hpp"
using morton::Morton;

Morton<2, 32> m = Morton<2, 32>::encode(x, y);  // interleave (PDEP)
m.inc(0);          // x += 1, in place, leaving y's interleaved bits untouched
m.add(1, 10);      // y += 10
auto n = m.neighbor(0, -1);   // the cell at (x-1, y+10)
auto [xx, yy] = m.decode();   // back to coordinates (PEXT)
```

## Why this is faster

Holding a Morton code and needing a neighbour is extremely common: stencils,
finite-difference sweeps, octree/quadtree traversal, ray marching, spatial
hashing. The conventional approach decodes to coordinates, changes one, and
re-encodes. This library mutates the relevant axis' bits in place.

Measured on this machine (Intel x86-64 with BMI2, g++ 14, `-O3`), 2D 32-bit:

| operation | this library | conventional (decode+encode) | speedup |
|---|--:|--:|--:|
| encode | 1398 Mops/s | 1400 Mops/s (libmorton) | parity |
| decode | 2040 Mops/s | 2012 Mops/s (libmorton) | parity |
| single neighbour (+1 on an axis) | **2721 Mops/s** | 1110 Mops/s | **2.5×** |
| 4-neighbour stencil (scattered) | **6718 Mops/s** | 2245 Mops/s | **3.0×** |

Encode/decode are at parity with [libmorton](https://github.com/Forceflow/libmorton)
(both use the hardware `PDEP`/`PEXT` instructions). The win is in the arithmetic
operations, which libmorton and similar libraries do not provide.

> **Honest caveat (from profiling):** for a *dense, contiguous* full-region
> sweep, simply re-encoding each `(x, y)` from scratch is actually faster than
> the arithmetic walk, because independent encodes pipeline while the arithmetic
> walk has a loop-carried dependency. The arithmetic advantage is for
> **scattered / data-dependent** neighbour access, where you already hold a code
> and cannot batch. See [`docs/EVALUATION.md`](docs/EVALUATION.md).

## Design

- **Header-only**, C++17, no dependencies. `#include "morton/morton.hpp"`.
- `Morton<Dim, Bits>` interleaves `Dim` coordinates of `Bits` bits each; the code
  is stored in the smallest unsigned integer that fits. `Dim * Bits ≤ 64` uses a
  built-in integer; **65–128 bits use `__uint128_t`** where available, so 3D
  32-bit (`Morton3D32`, 96-bit) and 2D 64-bit (`Morton2D64`, 128-bit) work too.
- Encode/decode use **BMI2 `PDEP`/`PEXT`** when compiled with `-mbmi2`, with a
  portable software fallback otherwise. Both paths are tested for agreement.
- **`constexpr`**: the software path runs at compile time (selected via
  `__builtin_is_constant_evaluated()`), so you can build lookup tables in
  `constexpr`.
- Arithmetic wraps modulo `2^Bits` per axis (branchless); `add_sat`/`sub_sat`/
  `try_add`/`try_sub` clamp or refuse instead of wrapping.

### Headers

| header | provides |
|---|---|
| `morton/morton.hpp` | `Morton<Dim,Bits>`: encode/decode, axis arithmetic, saturating ops, `face_neighbors`/`all_neighbors`, `ancestor`/`child` hierarchy, Z-order `++/--` |
| `morton/iterate.hpp` | `for_each_in_box` (row-major), `for_each_in_box_zorder` (Z-order), and the Tropf-Herzog range-search pair `bigmin_in_box` / `litmax_in_box` |
| `morton/batch.hpp` | vectorised bulk `add`/`sub`/`step`/`encode` over arrays (AVX2 auto-vectorised) |
| `morton/wide_uint.hpp` | fixed-width word-array unsigned backing codes wider than 128 bits (pulled in automatically) |

Convenience aliases: `Morton2D32`, `Morton2D16`, `Morton3D21`, `Morton3D16`,
`Morton3D32`, `Morton2D64`. Codes wider than 128 bits (up to `MORTON_MAX_BITS`,
default 256) work too, e.g. `Morton<3,64>` (192-bit), `Morton<2,128>` (256-bit).

A linear **octree/quadtree** built on this library lives in the sibling
[`octree/`](octree/) project (`morton_octree::Octree`) — it is being split into
its own repository; see [octree/PLAN.md](octree/PLAN.md).

### GPU (Kokkos)

A portable **Kokkos** backend lives in
[`include/morton/kokkos.hpp`](include/morton/kokkos.hpp) (`morton::kokkos`).
Because the core marks its functions with `MORTON_HD` (which defers to
`KOKKOS_FUNCTION`), the device kernels run the *same* `Morton<Dim,Bits>` code as
the CPU — no separate implementation — on **any** Kokkos backend (CUDA / HIP /
OpenMP / Serial); the backend and architecture come from the Kokkos build, not
the source. It offers `Kokkos::View`-based bulk ops (`encode2/3`, `decode2/3`,
per-axis `add`/`sub`/`step`) plus `*_host` raw-pointer convenience wrappers. On
an RTX 5080 (CUDA backend), 2D-32 encode hits **~51,000 Mops/s** when data is
resident on the GPU (~33× one CPU core); a one-shot call on host data is
PCIe-transfer-bound and no faster than the CPU, so the GPU pays off only when
codes live on-device across a pipeline. Opt in with `-DMORTON_ENABLE_KOKKOS=ON`
(see [DEVELOPMENT.md](DEVELOPMENT.md)).

## Build, test, benchmark

Full per-platform setup (incl. Python wheels and the GPU build) is in
[DEVELOPMENT.md](DEVELOPMENT.md). Quick start:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure   # or: ./build/tests/morton_tests
./build/benchmarks/morton_bench
```

As a header-only dependency in another CMake project — either vendored:

```cmake
add_subdirectory(morton_arithmetic)
target_link_libraries(your_target PRIVATE morton::morton)
```

…or installed and found via `find_package` (an `install` exports a CMake
package config; Conan recipe in `conanfile.py`, vcpkg port in
`packaging/vcpkg/morton/`):

```cmake
find_package(morton CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE morton::morton)
```

The test suite cross-checks every encode/decode against both a bit-by-bit
reference and libmorton, validates all arithmetic exhaustively on small grids
and against `decode/modify/encode` on large random inputs, and validates the
region iterators against brute force (>1M assertions).

## Python

A vectorised NumPy interface (pure `ctypes`, no pybind11/native build deps
beyond a C++ compiler), published on PyPI as **`peclet-morton`** and imported as
**`peclet.morton`** (part of the `peclet` suite namespace). Install as a wheel
with scikit-build-core:

```bash
pip install .            # builds the extension and bundles it into the wheel
python -m pytest bindings/python/tests -q
```

```python
import numpy as np, peclet.morton as ma
x = np.arange(1000, dtype=np.uint32); y = x * 2
codes = ma.encode(x, y, bits=32)
shifted = ma.shift(codes, axis=0, delta=+1, dims=2, bits=32)  # +1 in x, no decode
print(ma.decode(shifted, dims=2, bits=32))
```

For a quick dev loop without installing, drop the `.so` next to the package and
use `PYTHONPATH`:

```bash
cmake --build build --target mortonarith_c
PYTHONPATH=bindings/python python3 -m pytest bindings/python/tests -q
```

Every call runs over whole arrays in compiled code. The arithmetic `shift` is
~3× faster than decode+encode in Python too (see `bindings/python/bench_python.py`).

## Repository layout

```
include/morton/      the library: morton.hpp, iterate.hpp, batch.hpp, simd.hpp,
                     wide_uint.hpp, kokkos.hpp (portable GPU backend)
tests/               doctest suite (encode/decode, arithmetic, constexpr, wide,
                     neighbours, batch, iterate)
benchmarks/          C++ micro-benchmarks (vs libmorton) + batch/SIMD benchmark
bindings/python/     ctypes + NumPy wrapper and pytest tests
pyproject.toml       scikit-build-core wheel build + cibuildwheel config
conanfile.py         Conan recipe;  packaging/vcpkg/  vcpkg port
cmake/               CMake package-config template (find_package(morton))
.github/workflows/   ci.yml (build matrix) + release.yml (PyPI wheels on tag)
docs/                EVALUATION.md, ROADMAP.md, HILBERT_GPU_NOTES.md, Doxyfile
octree/              sibling project: linear octree on this library (being split out)
legacy/              the original arbitrary-width BitArray + octree prototype
third_party/         vendored doctest and libmorton (tests/benchmarks only)
```

The `legacy/` directory holds the project's original prototype: an
arbitrary-width `BitArray`, a wide-code `Morton`, and an octree built on it. The
new core supersedes it for codes that fit in 64 bits; see the roadmap for how
the octree could be ported.

## License

MIT — see [LICENSE](LICENSE).
