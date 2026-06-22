# Development setup

How to build, test, and develop the library on a fresh machine. The repository
is self-contained — clone it and follow the section you need. Nothing here is
specific to one workstation; ISA features are detected at configure time.

## Prerequisites

| Task | Needs |
|---|---|
| C++ core: build + test + benchmark | CMake ≥ 3.16, a C++17 compiler (GCC ≥ 9, Clang ≥ 9, or MSVC 2019+) |
| Python bindings / wheel | Python ≥ 3.8, `pip`; build pulls `scikit-build-core` + `numpy` |
| API docs | Doxygen |
| GPU (Kokkos) backend | a Kokkos install (CUDA / HIP / OpenMP) on `CMAKE_PREFIX_PATH` (see below) |

The core is **header-only and dependency-free**; `third_party/` vendors doctest
and libmorton for tests/benchmarks only.

## C++ core

```bash
git clone -b morton-arith-library \
  git@github.com:computational-chemical-engineering/morton_artithmetic.git
cd morton_artithmetic

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure        # or ./build/tests/morton_tests
./build/benchmarks/morton_bench                    # vs libmorton
./build/benchmarks/morton_bench_batch              # SIMD / batch
cmake --build build --target docs                  # Doxygen -> docs/doxygen/html
```

Useful options: `-DMORTON_ENABLE_BMI2=OFF` (force the portable software path),
`-DMORTON_BUILD_BENCHMARKS=OFF`, `-DMORTON_BUILD_BINDINGS=OFF`. A non-BMI2 build
must also pass — it exercises the software fallback:

```bash
cmake -S . -B build_nobmi2 -DMORTON_ENABLE_BMI2=OFF && cmake --build build_nobmi2 -j
```

### ISA / platform notes

- **BMI2 (`pdep`/`pext`)** is auto-detected (`check_cxx_compiler_flag`) and used
  for codes ≤ 64 bits. Fast on all modern Intel (Haswell 2013+) and AMD Zen3+;
  microcoded/slow on AMD Zen1/Zen2 — build with `-DMORTON_ENABLE_BMI2=OFF` there
  if you measure a regression.
- **`__int128`** (codes 65–128 bits) is used automatically on GCC/Clang x86-64;
  MSVC falls back to the `wide_uint` word-array backend. Codes > 128 bits always
  use `wide_uint` (cap `MORTON_MAX_BITS`, default 256; override with
  `-DMORTON_MAX_BITS=...` or `add_compile_definitions`).
- **AVX-512**: not required. `morton/batch.hpp` runtime-dispatches the bulk
  encode/decode/add to explicit AVX-512 kernels (`morton/simd.hpp`) when the CPU
  has AVX-512F, else the auto-vectorised scalar path runs. The kernels use a
  per-function `target("avx512f")` attribute, so no global `-mavx512f` is needed.

### Runtime dispatch (one portable binary)

`-DMORTON_ENABLE_RUNTIME_DISPATCH=ON` (pair with `-DMORTON_ENABLE_BMI2=OFF`)
builds a binary that contains **no** `-mbmi2` codegen globally but still uses
PDEP/PEXT at runtime when the CPU has BMI2 (software fallback otherwise), via
`target`-attribute helpers + a CPUID check. This is what the redistributable
wheels build, so a single wheel is SIGILL-safe on old CPUs yet BMI2/AVX-512-fast
on new ones. It self-disables under `-mbmi2`, on non-x86, and under MSVC/CUDA.

```bash
cmake -S . -B build_rt -DCMAKE_BUILD_TYPE=Release \
  -DMORTON_ENABLE_BMI2=OFF -DMORTON_ENABLE_RUNTIME_DISPATCH=ON
cmake --build build_rt -j
```

### Testing AVX-512 without AVX-512 hardware (Intel SDE)

Most dev machines and CI runners lack AVX-512, so validate the AVX-512 path
under [Intel SDE](https://www.intel.com/content/www/us/en/developer/articles/tool/software-development-emulator.html),
a userspace functional emulator (it runs on AMD too — no special hardware).
Download the Linux `.tar.xz`, extract it, and run the test binary under it:

```bash
# (one-time) fetch + unpack SDE somewhere, e.g. ~/sde
sde64 -skx -- ./build/tests/morton_tests        # Skylake-X: AVX-512 batch path runs
# the runtime-dispatch build across BMI2-present / -absent CPUs:
sde64 -snb -- ./build_rt/tests/morton_tests      # Sandy Bridge: no BMI2/AVX-512 -> software
sde64 -hsw -- ./build_rt/tests/morton_tests      # Haswell:      BMI2, no AVX-512 -> PDEP
sde64 -skx -- ./build_rt/tests/morton_tests      # Skylake-X:    BMI2 + AVX-512
```

The test suite compares the SIMD/dispatch output bit-for-bit against the scalar
reference, so a green run under `-skx` proves the AVX-512 kernels are correct.
CI does exactly this (the `cpp-avx512-sde` job), so AVX-512 is covered without a
GPU/AVX-512 runner. (Emulated *timing* is meaningless — benchmark on real
silicon; `morton_bench_batch` prints the active path and a bulk-encode compare.)

## Python bindings

```bash
pip install .                          # scikit-build-core builds + bundles the .so
python -m pytest bindings/python/tests -q
```

Dev loop without installing:

```bash
cmake --build build --target mortonarith_c     # drops .so into bindings/python/mortonarith/
PYTHONPATH=bindings/python python -m pytest bindings/python/tests -q
```

Redistributable wheels (`cibuildwheel`, see `pyproject.toml`) build the
**portable** software path (`MORTON_ENABLE_BMI2=OFF`) so they never SIGILL on a
CPU without BMI2. Source installs keep BMI2 on.

## Using the library as a dependency

```cmake
# vendored:
add_subdirectory(morton_arithmetic)
target_link_libraries(app PRIVATE morton::morton)

# installed (cmake --install build --prefix <p>):
find_package(morton CONFIG REQUIRED)
target_link_libraries(app PRIVATE morton::morton)
```

Conan recipe: `conanfile.py`. vcpkg port: `packaging/vcpkg/morton/`.

## GPU (Kokkos) backend — `include/morton/kokkos.hpp`

The portable Kokkos backend (`morton::kokkos`) reuses the core code path
(`MORTON_HD` → `KOKKOS_FUNCTION`), so there is nothing GPU-specific to maintain
in the math, and it runs on **any** Kokkos backend — CUDA, HIP, OpenMP or
Serial. The backend and target architecture come from the Kokkos build it is
linked against, not from this repo. (This replaces the retired raw-CUDA backend;
the last raw-CUDA tree is preserved under the `pre-cuda-retirement` git tag.)

Kokkos is consumed the same way as the rest of the `peclet` suite:
`find_package(Kokkos CONFIG)` against a prefix on `CMAKE_PREFIX_PATH`, provided
either by a cluster module (`module load Kokkos`) or by the suite bootstrap:

```bash
# one-time, from the suite root: build Kokkos for the backend you want
../tools/bootstrap_deps.sh host-openmp        # or nvidia-cuda / lumi-hip

# enable the morton Kokkos backend, point it at that prefix
cmake -S . -B build_kokkos -DMORTON_ENABLE_KOKKOS=ON \
      -DCMAKE_PREFIX_PATH="$PWD/../extern/install/host-openmp"
cmake --build build_kokkos -j
ctest --test-dir build_kokkos --output-on-failure   # morton_kokkos_tests
./build_kokkos/benchmarks/morton_bench_kokkos
```

For the CUDA backend, point `CMAKE_PREFIX_PATH` at `extern/install/nvidia-cuda`
and put `nvcc` on `PATH` (e.g. `export PATH=/usr/local/cuda-13.2/bin:$PATH`);
Kokkos routes the device sources through the launch compiler automatically — the
sources stay plain `.cpp`. The test cross-checks every device result against the
scalar `Morton<>` reference bit-for-bit on all four binding layouts
(`(2,32) (2,16) (3,21) (3,16)`).

GPU CI needs a GPU runner (self-hosted), so the CUDA/HIP backends are not part of
the default GitHub Actions matrix (the OpenMP backend is exercisable anywhere).

## Repository map

See the "Repository layout" section of [README.md](README.md). In short:
`include/morton/` is the library (incl. the `kokkos.hpp` GPU backend); `tests/`
(with `tests/kokkos/`), `benchmarks/`, `bindings/python/`, `octree/` (sibling
project), `docs/`, `packaging/`, `legacy/`.
