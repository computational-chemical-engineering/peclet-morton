# Evaluation: does this library contribute anything new?

Short answer: **the individual ideas are not new, but the packaging is useful and,
as far as I can tell, not offered as a unit by any popular library.** The honest
picture, including where the speed advantage does *and does not* materialise, is
below.

## 1. The landscape of existing Morton libraries

| Library | Lang | Encode/decode | Per-axis arithmetic | Z-order range query | Notes |
|---|---|---|---|---|---|
| [libmorton](https://github.com/Forceflow/libmorton) | C++ header-only | ✅ (LUT + BMI2 + AVX512) | ❌ | ❌ | The de-facto standard. Fastest encode/decode; nothing else. |
| [morton-nd](https://github.com/morton-nd/morton-nd) | C++17 header-only | ✅ (constexpr, N-dim, BMI2) | ❌ | ❌ | Elegant arbitrary-dimension encode/decode. No arithmetic. |
| [pymorton](https://github.com/trevorprater/pymorton) | Python | ✅ | ❌ (only via decode/encode) | ❌ | Pure Python, slow. |
| `numpy`-based snippets / [morton-py] | Python | ✅ (vectorised magic bits) | ❌ | ❌ | Encode/decode only. |
| **this library** | C++17 header-only + Python | ✅ (BMI2 + sw fallback) | ✅ inc/dec/add/sub/neighbour | ✅ BIGMIN row-major + Z-order | The combination is the contribution. |

So on encode/decode alone this library has **no advantage** over libmorton — it
is at parity (same `PDEP`/`PEXT` instructions), and libmorton additionally has an
AVX-512 batch path this library does not. If all you need is encode/decode, use
libmorton.

## 2. Is the arithmetic itself novel?

No. Doing integer arithmetic on a single dilated/interleaved coordinate while
leaving the others in place is **dilated integer arithmetic**, known since at
least the 1990s:

- D. S. Wise & J. Frens, *"Morton-order matrices…"* — dilated arithmetic for
  quadtree matrix layouts.
- L. Stocco & G. Schrack, *"Integer dilation and contraction for quadtrees and
  octrees"* (1995), and Schrack's *"Finding neighbors of equal size in linear
  quadtrees and octrees in constant time"* (1992) — this is exactly the
  `(code | ~mask) + 1` neighbour trick used here.
- H. Tropf & H. Herzog, *"Multidimensional Range Search in Dynamically Balanced
  Trees"* (1981) — the LITMAX/BIGMIN algorithm used for the Z-order range
  iterator.

The legacy prototype in `legacy/` already used the carry-propagation idea (on an
arbitrary-width bit array). What this library adds over that prototype is making
it **fast** (single-word, BMI2, branchless, O(1) instead of O(width) ripple loops)
and **general** (any `Dim`, any `Bits ≤ 64`, both 2D and 3D, software fallback).

## 3. So what *is* the contribution?

A single, modern, header-only, dependency-free C++17 library that combines, with
a tested and benchmarked implementation:

1. **Hardware-accelerated encode/decode** (parity with the best),
2. **O(1) per-axis arithmetic and neighbour finding** (which the popular
   encode/decode libraries simply do not expose), and
3. **Z-order range iteration** (BIGMIN), validated against brute force,

plus a **vectorised NumPy binding** that carries the same arithmetic advantage
into Python — where the gap is widest because the Python ecosystem's Morton
options are encode/decode-only and often pure-Python.

To my knowledge no widely-used library ships (1)+(2)+(3) together. That is a
real, if modest, gap to fill.

## 4. Where the speedup is real — and where it isn't

From `benchmarks/morton_bench` (x86-64, BMI2, g++ 14 `-O3`, 2D 32-bit):

```
encode (this library, PDEP)        1398 Mops/s   ┐ parity with libmorton
encode (libmorton)                 1400 Mops/s   ┘
decode (this library, PEXT)        2040 Mops/s   ┐ parity with libmorton
decode (libmorton)                 2012 Mops/s   ┘

neighbour +1 (arithmetic)          2721 Mops/s   ┐ 2.5×
neighbour +1 (decode+encode)       1110 Mops/s   ┘

4-neighbour stencil (arithmetic)   6718 Mops/s   ┐ 3.0×
4-neighbour stencil (decode+enc.)  2245 Mops/s   ┘

dense sweep (for_each_in_box)       895 Mops/s   ┐ arithmetic LOSES here
dense sweep (per-cell encode)      4510 Mops/s   ┘ (5×, the other way)
```

The decisive finding from profiling:

- **Scattered / data-dependent neighbour access wins big (2.5–3×).** You hold a
  code, you need an adjacent code, you cannot batch. This is the octree/stencil/
  spatial-hash access pattern, and it is where the library earns its keep.
- **Dense contiguous sweeps lose.** When you are enumerating an entire region,
  re-encoding each cell from scratch is *faster* than walking with arithmetic,
  because the independent `PDEP`s pipeline across the superscalar core while the
  arithmetic walk is a serial dependency chain (each cell depends on the previous
  one) with a per-cell carry branch. `for_each_in_box` is therefore offered for
  convenience and for the Z-order variant, not as the fast path for dense fills.

This nuance is the kind of thing that only shows up under measurement, and it
shapes the honest pitch: **this is a navigation library, not a faster encoder.**

## 5. Recommendation

The library is worth keeping and polishing as a **niche, well-scoped tool** for
code that *traverses* Z-ordered data (octrees/quadtrees, stencil solvers on
space-filling-curve layouts, spatial indexes). It should not be marketed as a
faster encoder — libmorton owns that — but as "libmorton + the arithmetic and
range-query operations it deliberately leaves out." See [`ROADMAP.md`](ROADMAP.md)
for the concrete plan.
