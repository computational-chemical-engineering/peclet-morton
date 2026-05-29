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
| GPU (CUDA) backend | an NVIDIA GPU + a CUDA toolkit (see below) |

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

## GPU (CUDA) backend — `cuda/`

Only needed if the machine has an NVIDIA GPU. The kernels reuse the core code
path (`__host__ __device__`), so there is nothing GPU-specific to maintain in
the math.

### Preferred: a normal CUDA Toolkit (`nvcc`)

Install the CUDA Toolkit (≥ 12.8 for Blackwell / sm_120; otherwise match your
GPU). Then:

```bash
cmake -S cuda -B cuda/build -DCMAKE_CUDA_ARCHITECTURES=90
cmake --build cuda/build -j
ctest --test-dir cuda/build --output-on-failure
./cuda/build/morton_cuda_bench
```

`CMAKE_CUDA_ARCHITECTURES=90` builds Hopper PTX, which the driver JITs onto
newer cards (e.g. sm_120 / RTX 50-series). Set it to your GPU's native arch for
a cubin (e.g. `120` with CUDA ≥ 12.8) to skip JIT.

### Fallback: clang as the CUDA compiler

Use this when you don't have a full `nvcc` (e.g. only the pip redistributable
pieces) or hit an `nvcc` vs. system-glibc header clash (the `cospi/sinpi`
"exception specification" errors with very new glibc + older CUDA). clang does
the CUDA front-end itself and needs only `ptxas` + `fatbinary` + `nvlink` +
`libdevice` + headers + `cudart`, arranged in a toolkit-shaped directory.

Assembling that directory from pip wheels + the NVIDIA redistributable archive:

```bash
python -m venv /tmp/cudaenv && /tmp/cudaenv/bin/pip install -U pip
/tmp/cudaenv/bin/pip install \
  nvidia-cuda-nvcc-cu12 nvidia-cuda-runtime-cu12 nvidia-curand-cu12 nvidia-cuda-cccl-cu12
NV=/tmp/cudaenv/lib/python*/site-packages/nvidia

# the nvcc pip wheel ships only ptxas; get nvcc/fatbinary/nvlink/cicc from the
# matching redistributable archive (version must match the nvcc wheel):
curl -sSL -o /tmp/nvcc.tar.xz \
  https://developer.download.nvidia.com/compute/cuda/redist/cuda_nvcc/linux-x86_64/cuda_nvcc-linux-x86_64-12.9.86-archive.tar.xz
mkdir -p /tmp/nvcc_pkg && tar -xJf /tmp/nvcc.tar.xz -C /tmp/nvcc_pkg --strip-components=1

# build a CUDA_PATH that clang understands
mkdir -p /tmp/cuda/bin /tmp/cuda/nvvm/libdevice /tmp/cuda/include /tmp/cuda/lib64
ln -sf /tmp/nvcc_pkg/bin/{ptxas,fatbinary,nvlink,cicc} /tmp/cuda/bin/
ln -sf $NV/cuda_nvcc/nvvm/libdevice/libdevice.10.bc /tmp/cuda/nvvm/libdevice/
cp -rs $NV/cuda_runtime/include/* $NV/cuda_nvcc/include/* $NV/curand/include/* \
       $NV/cuda_cccl/include/* /tmp/cuda/include/ 2>/dev/null
ln -sf $NV/cuda_runtime/lib/libcudart.so.12 /tmp/cuda/lib64/libcudart.so
ln -sf $NV/cuda_runtime/lib/libcudart.so.12 /tmp/cuda/lib64/libcudart.so.12

# build + run (sm_90 PTX JITs onto newer GPUs)
CUDA_PATH=/tmp/cuda GPU_ARCH=sm_90 ./cuda/build_clang.sh
```

`cuda/build_clang.sh` compiles and runs the tests and benchmark; it only needs
`CUDA_PATH` (and optionally `GPU_ARCH`, `CXX`). Run binaries with
`LD_LIBRARY_PATH=$CUDA_PATH/lib64`.

> The exact `cuda_nvcc` archive name (`...-12.9.86-...`) must match the installed
> `nvidia-cuda-nvcc-cu12` version; list options in
> `https://developer.download.nvidia.com/compute/cuda/redist/redistrib_<ver>.json`.

GPU CI needs a GPU runner (self-hosted), so the CUDA backend is not part of the
default GitHub Actions matrix.

## Repository map

See the "Repository layout" section of [README.md](README.md). In short:
`include/morton/` is the library; `tests/`, `benchmarks/`, `bindings/python/`,
`cuda/`, `octree/` (sibling project), `docs/`, `packaging/`, `legacy/`.
