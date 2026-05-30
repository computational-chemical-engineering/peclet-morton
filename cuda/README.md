# morton-cuda

CUDA backend for the [morton-arithmetic](../README.md) library.

The per-element work is the ordinary `Morton<Dim,Bits>` code path. The core
header marks its functions `__host__ __device__` (via `MORTON_HD`) when compiled
with a CUDA/HIP compiler, so the GPU kernels run the **exact, already-tested
logic of the CPU library** — there is no second implementation to keep in sync.
Encode/decode on the GPU use the software (PDEP-free) bit path, since there is
no PDEP/PEXT on NVIDIA hardware.

## API (`morton_cuda/morton_cuda.cuh`, namespace `morton::cuda`)

- `encode2_device` / `encode3_device` / `decode2_device` / `add_device` — operate
  on **device pointers** the caller owns (compose into pipelines, choose streams).
- `encode2_host` / `decode2_host` / `encode3_host` / `add_host` — convenience
  wrappers that allocate, copy host↔device, launch and free (one-shot calls).

```cpp
#include "morton_cuda/morton_cuda.cuh"
// x, y, out are host arrays of length n
morton::cuda::encode2_host<32>(x, y, out, n);            // -> Morton codes
morton::cuda::add_host<2,32>(out, out, n, /*axis*/0, +1); // move +1 in x on the GPU
```

## Build

With a normal CUDA Toolkit (`nvcc`) via CMake:

```bash
cmake -S cuda -B cuda/build -DCMAKE_CUDA_ARCHITECTURES=90
cmake --build cuda/build -j
ctest --test-dir cuda/build --output-on-failure
```

Use your GPU's architecture; `90` (Hopper) PTX JITs onto newer cards such as
sm_120 / RTX 50-series via the driver.

Without a full `nvcc` (e.g. only ptxas + libdevice + cudart from the pip
redistributables, or to dodge nvcc/glibc header clashes), build with **clang as
the CUDA compiler**:

```bash
CUDA_PATH=/path/to/cuda GPU_ARCH=sm_90 ./cuda/build_clang.sh
```

## Performance (RTX 5080, 2D 32-bit encode, ~34M points)

```
GPU device-only        ~51000 Mops/s   (data resident on the GPU)
GPU host round-trip     ~1200 Mops/s   (incl. H2D + D2H copies)
CPU batch (1 core)      ~1550 Mops/s
```

The device-only number (~33× one CPU core) is the win when codes **live on the
GPU** across a pipeline. A one-shot call on host data is **transfer-bound** —
the PCIe copies dominate and it is no faster than the CPU. So put the GPU at the
*start* of a chain (encode → sort → query → … all on device), not as a drop-in
for a single host call. This mirrors the CPU finding that bulk ops are
memory-bound; here the bottleneck just moves to the bus.

## Status / future

Implemented: encode/decode (2D/3D) and per-axis arithmetic. Natural next steps
(see [../docs/HILBERT_GPU_NOTES.md](../docs/HILBERT_GPU_NOTES.md)): a Z-order
**radix sort** (the usual reason to be on the GPU), a NumPy-device Python entry
point, and SYCL/HIP portability. GPU CI needs a GPU runner.
