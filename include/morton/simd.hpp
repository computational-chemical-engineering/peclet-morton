// morton/simd.hpp
//
// Explicit AVX-512 batch kernels for the bulk encode/decode/arithmetic paths.
// Unlike the auto-vectorised loops in batch.hpp, these are hand-written with
// AVX-512 intrinsics ("magic-bits" bit spreading, libmorton-style) and process
// eight 64-bit codes per iteration.
//
// They are compiled with a per-function `__attribute__((target("avx512f")))`,
// so the translation unit does *not* need -mavx512f and the rest of the binary
// stays portable. batch.hpp dispatches to them at runtime only when
// detail::cpu_has_avx512f() is true, so a single binary uses AVX-512 on a
// capable CPU and the scalar/AVX2 path elsewhere. Codes that are not 64 bits
// wide keep the scalar path.
//
// SPDX-License-Identifier: MIT

#ifndef MORTON_SIMD_HPP
#define MORTON_SIMD_HPP

#include <cstddef>
#include <cstdint>

#include "morton/morton.hpp"

#if MORTON_X86

#include <immintrin.h>

#define MORTON_AVX512_TARGET __attribute__((target("avx512f")))

namespace morton {
namespace detail {
namespace avx512 {

// ---- 2D (one zero-bit gap): spread 32-bit -> 64-bit, and its inverse --------

MORTON_AVX512_TARGET inline __m512i spread2(__m512i x) {
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_slli_epi64(x, 16)),
                         _mm512_set1_epi64(0x0000FFFF0000FFFFLL));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_slli_epi64(x, 8)),
                         _mm512_set1_epi64(0x00FF00FF00FF00FFLL));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_slli_epi64(x, 4)),
                         _mm512_set1_epi64(0x0F0F0F0F0F0F0F0FLL));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_slli_epi64(x, 2)),
                         _mm512_set1_epi64(0x3333333333333333LL));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_slli_epi64(x, 1)),
                         _mm512_set1_epi64(0x5555555555555555LL));
    return x;
}

MORTON_AVX512_TARGET inline __m512i compact2(__m512i x) {
    x = _mm512_and_si512(x, _mm512_set1_epi64(0x5555555555555555LL));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi64(x, 1)),
                         _mm512_set1_epi64(0x3333333333333333LL));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi64(x, 2)),
                         _mm512_set1_epi64(0x0F0F0F0F0F0F0F0FLL));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi64(x, 4)),
                         _mm512_set1_epi64(0x00FF00FF00FF00FFLL));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi64(x, 8)),
                         _mm512_set1_epi64(0x0000FFFF0000FFFFLL));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi64(x, 16)),
                         _mm512_set1_epi64(0x00000000FFFFFFFFLL));
    return x;
}

// ---- 3D (two zero-bit gaps): spread 21-bit -> 63-bit, and its inverse -------

MORTON_AVX512_TARGET inline __m512i spread3(__m512i x) {
    x = _mm512_and_si512(x, _mm512_set1_epi64(0x1FFFFFLL));  // keep low 21 bits
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_slli_epi64(x, 32)),
                         _mm512_set1_epi64(0x1F00000000FFFFLL));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_slli_epi64(x, 16)),
                         _mm512_set1_epi64(0x1F0000FF0000FFLL));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_slli_epi64(x, 8)),
                         _mm512_set1_epi64(0x100F00F00F00F00FLL));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_slli_epi64(x, 4)),
                         _mm512_set1_epi64(0x10C30C30C30C30C3LL));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_slli_epi64(x, 2)),
                         _mm512_set1_epi64(0x1249249249249249LL));
    return x;
}

MORTON_AVX512_TARGET inline __m512i compact3(__m512i x) {
    x = _mm512_and_si512(x, _mm512_set1_epi64(0x1249249249249249LL));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi64(x, 2)),
                         _mm512_set1_epi64(0x10C30C30C30C30C3LL));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi64(x, 4)),
                         _mm512_set1_epi64(0x100F00F00F00F00FLL));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi64(x, 8)),
                         _mm512_set1_epi64(0x1F0000FF0000FFLL));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi64(x, 16)),
                         _mm512_set1_epi64(0x1F00000000FFFFLL));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi64(x, 32)),
                         _mm512_set1_epi64(0x1FFFFFLL));
    return x;
}

// ---- encode / decode (2,32) and (3,21) --------------------------------------
// The scalar tail (n not a multiple of 8) reuses the core Morton<> path, so the
// results are bit-for-bit identical to the scalar functions.

MORTON_AVX512_TARGET inline void encode2(const std::uint32_t* x, const std::uint32_t* y,
                                         std::uint64_t* out, std::size_t n) {
    std::size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m512i vx = _mm512_cvtepu32_epi64(_mm256_loadu_si256((const __m256i*)(x + i)));
        __m512i vy = _mm512_cvtepu32_epi64(_mm256_loadu_si256((const __m256i*)(y + i)));
        __m512i code = _mm512_or_si512(spread2(vx), _mm512_slli_epi64(spread2(vy), 1));
        _mm512_storeu_si512((__m512i*)(out + i), code);
    }
    for (; i < n; ++i) out[i] = Morton<2, 32>::encode(x[i], y[i]).code();
}

MORTON_AVX512_TARGET inline void decode2(const std::uint64_t* in, std::uint32_t* x,
                                         std::uint32_t* y, std::size_t n) {
    std::size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m512i c = _mm512_loadu_si512((const __m512i*)(in + i));
        _mm256_storeu_si256((__m256i*)(x + i), _mm512_cvtepi64_epi32(compact2(c)));
        _mm256_storeu_si256((__m256i*)(y + i),
                            _mm512_cvtepi64_epi32(compact2(_mm512_srli_epi64(c, 1))));
    }
    for (; i < n; ++i) {
        auto a = Morton<2, 32>::from_code(in[i]).decode();
        x[i] = a[0];
        y[i] = a[1];
    }
}

MORTON_AVX512_TARGET inline void encode3(const std::uint32_t* x, const std::uint32_t* y,
                                         const std::uint32_t* z, std::uint64_t* out,
                                         std::size_t n) {
    std::size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m512i vx = _mm512_cvtepu32_epi64(_mm256_loadu_si256((const __m256i*)(x + i)));
        __m512i vy = _mm512_cvtepu32_epi64(_mm256_loadu_si256((const __m256i*)(y + i)));
        __m512i vz = _mm512_cvtepu32_epi64(_mm256_loadu_si256((const __m256i*)(z + i)));
        __m512i code = _mm512_or_si512(
            _mm512_or_si512(spread3(vx), _mm512_slli_epi64(spread3(vy), 1)),
            _mm512_slli_epi64(spread3(vz), 2));
        _mm512_storeu_si512((__m512i*)(out + i), code);
    }
    for (; i < n; ++i) out[i] = Morton<3, 21>::encode(x[i], y[i], z[i]).code();
}

MORTON_AVX512_TARGET inline void decode3(const std::uint64_t* in, std::uint32_t* x,
                                         std::uint32_t* y, std::uint32_t* z, std::size_t n) {
    std::size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m512i c = _mm512_loadu_si512((const __m512i*)(in + i));
        _mm256_storeu_si256((__m256i*)(x + i), _mm512_cvtepi64_epi32(compact3(c)));
        _mm256_storeu_si256((__m256i*)(y + i),
                            _mm512_cvtepi64_epi32(compact3(_mm512_srli_epi64(c, 1))));
        _mm256_storeu_si256((__m256i*)(z + i),
                            _mm512_cvtepi64_epi32(compact3(_mm512_srli_epi64(c, 2))));
    }
    for (; i < n; ++i) {
        auto a = Morton<3, 21>::from_code(in[i]).decode();
        x[i] = a[0];
        y[i] = a[1];
        z[i] = a[2];
    }
}

// ---- per-axis add / sub on 64-bit codes (any 64-bit-code layout) ------------
// Mirrors the masked-add identity in batch.hpp; the mask/increment are
// loop-invariant scalars supplied by the caller.

MORTON_AVX512_TARGET inline void add64(const std::uint64_t* in, std::uint64_t* out,
                                       std::size_t n, std::uint64_t notM, std::uint64_t dk,
                                       std::uint64_t mask, std::uint64_t keep) {
    const __m512i vnotM = _mm512_set1_epi64((long long)notM);
    const __m512i vdk = _mm512_set1_epi64((long long)dk);
    const __m512i vmask = _mm512_set1_epi64((long long)mask);
    const __m512i vkeep = _mm512_set1_epi64((long long)keep);
    std::size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m512i c = _mm512_loadu_si512((const __m512i*)(in + i));
        __m512i t = _mm512_add_epi64(_mm512_or_si512(c, vnotM), vdk);
        __m512i r = _mm512_or_si512(_mm512_and_si512(t, vmask), _mm512_and_si512(c, vkeep));
        _mm512_storeu_si512((__m512i*)(out + i), r);
    }
    for (; i < n; ++i) out[i] = (((in[i] | notM) + dk) & mask) | (in[i] & keep);
}

MORTON_AVX512_TARGET inline void sub64(const std::uint64_t* in, std::uint64_t* out,
                                       std::size_t n, std::uint64_t dk, std::uint64_t mask,
                                       std::uint64_t keep) {
    const __m512i vdk = _mm512_set1_epi64((long long)dk);
    const __m512i vmask = _mm512_set1_epi64((long long)mask);
    const __m512i vkeep = _mm512_set1_epi64((long long)keep);
    std::size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m512i c = _mm512_loadu_si512((const __m512i*)(in + i));
        __m512i t = _mm512_sub_epi64(_mm512_and_si512(c, vmask), vdk);
        __m512i r = _mm512_or_si512(_mm512_and_si512(t, vmask), _mm512_and_si512(c, vkeep));
        _mm512_storeu_si512((__m512i*)(out + i), r);
    }
    for (; i < n; ++i) out[i] = (((in[i] & mask) - dk) & mask) | (in[i] & keep);
}

}  // namespace avx512
}  // namespace detail
}  // namespace morton

#undef MORTON_AVX512_TARGET

#endif  // MORTON_X86

#endif  // MORTON_SIMD_HPP
