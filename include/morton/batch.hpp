// morton/batch.hpp
//
// Bulk array operations. The per-axis arithmetic is a short sequence of masked
// integer ops with a loop-invariant mask and increment, so these loops
// auto-vectorise (AVX2/AVX-512) under -O3: the compiler turns them into packed
// vpor/vpaddq/vpand. This is what gives the Python bindings near-native
// throughput, and what a SIMD-batch API would expose to C++ callers.
//
// SPDX-License-Identifier: MIT

#ifndef MORTON_BATCH_HPP
#define MORTON_BATCH_HPP

#include <cstddef>
#include <cstdint>

#include "morton/morton.hpp"

namespace morton {
namespace batch {

/// Add the same `k` to `axis` of every code (wraps). Vectorised.
/// `in` and `out` may alias.
template <unsigned Dim, unsigned Bits>
void add(const typename Morton<Dim, Bits>::code_type* in,
         typename Morton<Dim, Bits>::code_type* out, std::size_t n, unsigned axis,
         typename Morton<Dim, Bits>::coord_type k) {
    using M = Morton<Dim, Bits>;
    using code = typename M::code_type;
    const code mask = M::axis_mask(axis);
    const code notM = code(~mask);
    const code keep = M::field_mask & notM;
    const code dk = M::deposit(k, axis);  // loop-invariant
    for (std::size_t i = 0; i < n; ++i) {
        code c = in[i];
        out[i] = code((c | notM) + dk) & mask;
        out[i] |= c & keep;
    }
}

/// Subtract the same `k` from `axis` of every code (wraps). Vectorised.
template <unsigned Dim, unsigned Bits>
void sub(const typename Morton<Dim, Bits>::code_type* in,
         typename Morton<Dim, Bits>::code_type* out, std::size_t n, unsigned axis,
         typename Morton<Dim, Bits>::coord_type k) {
    using M = Morton<Dim, Bits>;
    using code = typename M::code_type;
    const code mask = M::axis_mask(axis);
    const code keep = M::field_mask & code(~mask);
    const code dk = M::deposit(k, axis);
    for (std::size_t i = 0; i < n; ++i) {
        code c = in[i];
        out[i] = (code((c & mask) - dk) & mask) | (c & keep);
    }
}

/// Step every code one cell along ±`axis` (dir>=0 -> +1, else -1). Vectorised.
template <unsigned Dim, unsigned Bits>
void step(const typename Morton<Dim, Bits>::code_type* in,
          typename Morton<Dim, Bits>::code_type* out, std::size_t n, unsigned axis, int dir) {
    if (dir >= 0)
        add<Dim, Bits>(in, out, n, axis, 1);
    else
        sub<Dim, Bits>(in, out, n, axis, 1);
}

/// Encode arrays of coordinates (2D) into codes.
template <unsigned Bits>
void encode2(const typename Morton<2, Bits>::coord_type* x,
             const typename Morton<2, Bits>::coord_type* y,
             typename Morton<2, Bits>::code_type* out, std::size_t n) {
    using M = Morton<2, Bits>;
    for (std::size_t i = 0; i < n; ++i) out[i] = M::encode(x[i], y[i]).code();
}

/// Encode arrays of coordinates (3D) into codes.
template <unsigned Bits>
void encode3(const typename Morton<3, Bits>::coord_type* x,
             const typename Morton<3, Bits>::coord_type* y,
             const typename Morton<3, Bits>::coord_type* z,
             typename Morton<3, Bits>::code_type* out, std::size_t n) {
    using M = Morton<3, Bits>;
    for (std::size_t i = 0; i < n; ++i) out[i] = M::encode(x[i], y[i], z[i]).code();
}

}  // namespace batch
}  // namespace morton

#endif  // MORTON_BATCH_HPP
