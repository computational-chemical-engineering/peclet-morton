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
- `Morton<Dim, Bits>` interleaves `Dim` coordinates of `Bits` bits each
  (`Dim * Bits ≤ 64`); the code is stored in the smallest unsigned integer that
  fits.
- Encode/decode use **BMI2 `PDEP`/`PEXT`** when compiled with `-mbmi2`, with a
  portable constexpr software fallback otherwise. Both paths are tested for
  agreement.
- Arithmetic wraps modulo `2^Bits` per axis (well-defined, branchless).
- `morton/iterate.hpp` adds region traversal: `for_each_in_box` (row-major,
  arithmetic) and `for_each_in_box_zorder` (Z-order, via the Tropf-Herzog
  BIGMIN range-search algorithm).

Convenience aliases: `Morton2D32`, `Morton2D16`, `Morton3D21`, `Morton3D16`.

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
beyond a C++ compiler) lives in `bindings/python`:

```bash
cmake --build build --target mortonarith_c        # builds the .so into the package
PYTHONPATH=bindings/python python3 -c "
import numpy as np, mortonarith as ma
x = np.arange(1000, dtype=np.uint32); y = x * 2
codes = ma.encode(x, y, bits=32)
shifted = ma.shift(codes, axis=0, delta=+1, dims=2, bits=32)  # +1 in x, no decode
print(ma.decode(shifted, dims=2, bits=32))
"
PYTHONPATH=bindings/python python3 -m pytest bindings/python/tests -q
```

Every call runs over whole arrays in compiled code. The arithmetic `shift` is
~3× faster than decode+encode in Python too (see `bindings/python/bench_python.py`).

## Repository layout

```
include/morton/      the library (morton.hpp, iterate.hpp)
tests/               doctest suite
benchmarks/          C++ micro-benchmarks (vs libmorton)
bindings/python/     ctypes + NumPy wrapper and pytest tests
docs/                EVALUATION.md (vs prior art) and ROADMAP.md
legacy/              the original arbitrary-width BitArray + octree prototype
third_party/         vendored doctest and libmorton (tests/benchmarks only)
```

The `legacy/` directory holds the project's original prototype: an
arbitrary-width `BitArray`, a wide-code `Morton`, and an octree built on it. The
new core supersedes it for codes that fit in 64 bits; see the roadmap for how
the octree could be ported.

## License

MIT — see [LICENSE](LICENSE).
