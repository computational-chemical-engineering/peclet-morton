# morton-arithmetic

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
| `morton/octree.hpp` | `Octree<Dim,Bits,T>`: linear octree/quadtree with point location, face neighbours and refinement built on the arithmetic core |

Convenience aliases: `Morton2D32`, `Morton2D16`, `Morton3D21`, `Morton3D16`,
`Morton3D32`, `Morton2D64`.

## Build, test, benchmark

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure   # or: ./build/tests/morton_tests
./build/benchmarks/morton_bench
```

As a header-only dependency in another CMake project:

```cmake
add_subdirectory(morton_arithmetic)
target_link_libraries(your_target PRIVATE morton::morton)
```

The test suite cross-checks every encode/decode against both a bit-by-bit
reference and libmorton, validates all arithmetic exhaustively on small grids
and against `decode/modify/encode` on large random inputs, and validates the
region iterators against brute force (>1M assertions).

## Python

A vectorised NumPy interface (pure `ctypes`, no pybind11/native build deps
beyond a C++ compiler). Install as a wheel with scikit-build-core:

```bash
pip install .            # builds the extension and bundles it into the wheel
python -m pytest bindings/python/tests -q
```

```python
import numpy as np, mortonarith as ma
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
include/morton/      the library: morton.hpp, iterate.hpp, batch.hpp, octree.hpp
tests/               doctest suite (encode/decode, arithmetic, constexpr, wide,
                     neighbours, octree, batch, iterate)
benchmarks/          C++ micro-benchmarks (vs libmorton) + batch/SIMD benchmark
bindings/python/     ctypes + NumPy wrapper and pytest tests
pyproject.toml       scikit-build-core wheel build
.github/workflows/   CI matrix (gcc/clang/MSVC x BMI2 x Debug/Release, Python, docs)
docs/                EVALUATION.md (vs prior art), ROADMAP.md, Doxyfile
legacy/              the original arbitrary-width BitArray + octree prototype
third_party/         vendored doctest and libmorton (tests/benchmarks only)
```

The `legacy/` directory holds the project's original prototype: an
arbitrary-width `BitArray`, a wide-code `Morton`, and an octree built on it. The
new core supersedes it for codes that fit in 64 bits; see the roadmap for how
the octree could be ported.

## License

MIT — see [LICENSE](LICENSE).
