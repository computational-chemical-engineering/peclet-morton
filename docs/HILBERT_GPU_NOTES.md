# Design notes: Hilbert curve and GPU batch kernels

Two larger, separable future efforts. These are notes, not commitments —
captured so the work can start without re-deriving the design.

## Hilbert curve support

### Why
The Hilbert curve has strictly better spatial locality than Morton/Z-order:
consecutive indices are always face-adjacent (no "long jumps" across the
domain), which improves cache/IO behaviour for range queries and tiling. Many
spatial-index and HPC workloads prefer it. It is the natural second curve for a
"space-filling-curve arithmetic" library.

### What carries over from the Morton core
- **Storage / width**: the `uint_for`/`wide_uint` machinery and the
  `Dim*Bits ≤ MORTON_MAX_BITS` framing are curve-agnostic — reuse as-is.
- **Range search**: BIGMIN/LITMAX have Hilbert analogues; the iterator structure
  in `iterate.hpp` (probe → skip gap → resume) is the same shape.
- **Comparisons / successor**: a Hilbert index is still a plain integer, so
  Z-order `++`/`--`, `operator<`, and the `std::map`-keying pattern transfer.

### What is genuinely different
- **encode/decode**: Hilbert needs the per-level rotation/reflection state
  machine (Skilling's transpose↔index algorithm is the cleanest: operate on the
  "transpose" form, then Gray-code in/out). No `PDEP`/`PEXT` shortcut — it is
  inherently sequential over the bit levels, so expect it to be slower than
  Morton encode and *not* at parity with libmorton-style table encoders.
- **Axis arithmetic — the headline feature does NOT transfer directly.** "Add 1
  to x" is O(1) on a Morton code because each axis owns a fixed set of bits.
  Hilbert indices interleave *and* rotate, so a single-axis step is not a fixed
  masked add. The honest options:
  1. decode → step coordinate → encode (the conventional path; no speedup), or
  2. precompute small per-(state,direction) transition tables to step the
     Hilbert index a level at a time — worth prototyping and benchmarking, but
     it will not be the "few-instructions" win Morton enjoys.

### Suggested shape
A separate `hilbert/` header (or a `Curve` policy template parameterising
`encode`/`decode`/`successor` while sharing storage and range-search code).
Keep Morton's arithmetic guarantees explicit so users don't assume Hilbert has
the same O(1) neighbour step. **Recommendation:** ship Hilbert as encode/decode
+ range-search first; treat fast Hilbert neighbour-stepping as a research item.

## GPU / SYCL batch kernels

> **Update: a portable Kokkos backend is now implemented** in
> `include/morton/kokkos.hpp` (`morton::kokkos`; it superseded the original
> raw-CUDA backend, retired at the `pre-cuda-retirement` tag), runs on CUDA / HIP
> / OpenMP, and follows the design below: the core is device-callable so kernels reuse
> the CPU code path, encode/decode use the software bit path, and the host API is
> array-shaped. Measured behaviour confirmed the prediction — device-resident
> throughput is ~33× one CPU core, but one-shot host calls are PCIe-transfer-
> bound. What remains from this section: the Z-order **radix sort**, the
> device-array Python entry point, and SYCL/HIP portability.

### Why
The bulk array ops in `batch.hpp` are embarrassingly parallel; large point sets
(graphics, particle sims, databases) encode/sort/query millions of codes.

### Design
- The arithmetic is pure integer ops on the code word(s) — directly portable to
  CUDA/SYCL/HIP. For `≤64`-bit codes a thread handles one code; `__uint128_t`
  maps to a 2×`u64` thread-local; `wide_uint` maps to a fixed `u64[W]`.
- `PDEP`/`PEXT` have no GPU equivalent, so GPU encode/decode uses the **software
  magic-bit/`spread_sw` path** (already the constexpr fallback — share it).
- Provide kernels for `encode`, `decode`, `add/step` (per-axis), and a sort
  (Z-order sort = radix sort on the code) which is the usual reason to go to GPU.
- Keep the host API NumPy-shaped (arrays in, arrays out) so the existing
  `bindings/` story extends: a `mortonarith.cuda` submodule mirroring the CPU
  functions, dispatching on array device.

### Practicalities
- Gate behind a build option; do not add a hard CUDA/SYCL dependency to the
  header-only core.
- Benchmark against the CPU `batch.hpp` path — for many workloads the bottleneck
  is the host↔device transfer, not the arithmetic (mirrors the CPU finding that
  bulk ops are already memory-bound). Measure before claiming a win.
- Needs CI with a GPU runner (or at least a SYCL CPU backend) before shipping;
  like AVX-512, untested SIMD/GPU code should not land.

## Relationship to the standalone decision

Both are candidates for **separate optional packages** depending on the morton
core (as the octree now is), keeping the core small, dependency-free, and easy
to audit.
