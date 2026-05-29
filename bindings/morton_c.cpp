// C ABI shim exposing the morton-arithmetic library to other languages
// (used by the ctypes-based Python package). All functions operate on whole
// arrays so the per-element Python overhead is amortised away and NumPy gets
// near-native throughput.
//
// The bulk encode/decode/add/sub go through morton::batch, which dispatches to
// AVX-512 (and runtime BMI2) on a capable CPU and the auto-vectorised scalar
// loop otherwise -- so a single portable .so adapts to the host. The C ABI
// always passes codes as uint64; configs whose code_type is uint64 (every one
// except 2D-16) take the batch path, the rest use the per-element core path.
//
// SPDX-License-Identifier: MIT

#include <array>
#include <cstddef>
#include <cstdint>
#include <type_traits>

#include "morton/batch.hpp"
#include "morton/iterate.hpp"
#include "morton/morton.hpp"

using morton::Morton;

namespace {

// Dispatch helpers as *templates* so `if constexpr` discards the branch that
// would not type-check for a given config (the uint64-ABI vs uint32-code_type
// mismatch on the 2D-16 layout). The discarded branch is never instantiated.

template <unsigned Dim, unsigned Bits>
constexpr bool code_is_u64 =
    std::is_same_v<typename Morton<Dim, Bits>::code_type, std::uint64_t>;

template <unsigned Bits, typename COORD>
inline void enc2(const COORD* x, const COORD* y, std::uint64_t* out, std::size_t n) {
    using M = Morton<2, Bits>;
    if constexpr (code_is_u64<2, Bits>)
        morton::batch::encode2<Bits>(x, y, out, n);
    else
        for (std::size_t i = 0; i < n; ++i) out[i] = M::encode(x[i], y[i]).code();
}

template <unsigned Bits, typename COORD>
inline void dec2(const std::uint64_t* code, COORD* x, COORD* y, std::size_t n) {
    using M = Morton<2, Bits>;
    if constexpr (code_is_u64<2, Bits>) {
        morton::batch::decode2<Bits>(code, x, y, n);
    } else {
        for (std::size_t i = 0; i < n; ++i) {
            M m = M::from_code(code[i]);
            x[i] = COORD(m.get(0));
            y[i] = COORD(m.get(1));
        }
    }
}

template <unsigned Bits>
inline void add2(const std::uint64_t* code, std::uint64_t* out, std::size_t n, unsigned axis,
                 std::int64_t k) {
    using M = Morton<2, Bits>;
    using C = typename M::coord_type;
    if constexpr (code_is_u64<2, Bits>) {
        if (k >= 0)
            morton::batch::add<2, Bits>(code, out, n, axis, C(std::uint64_t(k)));
        else
            morton::batch::sub<2, Bits>(code, out, n, axis, C(std::uint64_t(-k)));
    } else {
        for (std::size_t i = 0; i < n; ++i) {
            M m = M::from_code(code[i]);
            if (k >= 0)
                m.add(axis, C(std::uint64_t(k)));
            else
                m.sub(axis, C(std::uint64_t(-k)));
            out[i] = m.code();
        }
    }
}

template <unsigned Bits, typename COORD>
inline void enc3(const COORD* x, const COORD* y, const COORD* z, std::uint64_t* out,
                 std::size_t n) {
    using M = Morton<3, Bits>;
    if constexpr (code_is_u64<3, Bits>)
        morton::batch::encode3<Bits>(x, y, z, out, n);
    else
        for (std::size_t i = 0; i < n; ++i) out[i] = M::encode(x[i], y[i], z[i]).code();
}

template <unsigned Bits, typename COORD>
inline void dec3(const std::uint64_t* code, COORD* x, COORD* y, COORD* z, std::size_t n) {
    using M = Morton<3, Bits>;
    if constexpr (code_is_u64<3, Bits>) {
        morton::batch::decode3<Bits>(code, x, y, z, n);
    } else {
        for (std::size_t i = 0; i < n; ++i) {
            M m = M::from_code(code[i]);
            x[i] = COORD(m.get(0));
            y[i] = COORD(m.get(1));
            z[i] = COORD(m.get(2));
        }
    }
}

template <unsigned Bits>
inline void add3(const std::uint64_t* code, std::uint64_t* out, std::size_t n, unsigned axis,
                 std::int64_t k) {
    using M = Morton<3, Bits>;
    using C = typename M::coord_type;
    if constexpr (code_is_u64<3, Bits>) {
        if (k >= 0)
            morton::batch::add<3, Bits>(code, out, n, axis, C(std::uint64_t(k)));
        else
            morton::batch::sub<3, Bits>(code, out, n, axis, C(std::uint64_t(-k)));
    } else {
        for (std::size_t i = 0; i < n; ++i) {
            M m = M::from_code(code[i]);
            if (k >= 0)
                m.add(axis, C(std::uint64_t(k)));
            else
                m.sub(axis, C(std::uint64_t(-k)));
            out[i] = m.code();
        }
    }
}

}  // namespace

#define API extern "C" __attribute__((visibility("default")))

// ---------------------------------------------------------------------------
// 2D configurations
// ---------------------------------------------------------------------------
#define DEFINE_2D(SUFFIX, BITS, COORD)                                              \
    API void mortonarith_encode2_##SUFFIX(const COORD* x, const COORD* y,           \
                                          std::uint64_t* out, std::size_t n) {      \
        enc2<BITS>(x, y, out, n);                                                   \
    }                                                                               \
    API void mortonarith_decode2_##SUFFIX(const std::uint64_t* code, COORD* x,      \
                                          COORD* y, std::size_t n) {                \
        dec2<BITS>(code, x, y, n);                                                  \
    }                                                                               \
    API void mortonarith_add2_##SUFFIX(const std::uint64_t* code, std::uint64_t* out, \
                                       std::size_t n, unsigned axis,                 \
                                       std::int64_t k) {                             \
        add2<BITS>(code, out, n, axis, k);                                          \
    }                                                                               \
    API std::uint64_t mortonarith_box_count2_##SUFFIX(const COORD* lo,              \
                                                      const COORD* hi) {            \
        std::uint64_t c = 1;                                                        \
        for (int d = 0; d < 2; ++d) {                                               \
            if (hi[d] < lo[d]) return 0;                                            \
            c *= std::uint64_t(hi[d] - lo[d]) + 1;                                  \
        }                                                                           \
        return c;                                                                   \
    }                                                                               \
    API void mortonarith_box_zorder2_##SUFFIX(const COORD* lo, const COORD* hi,     \
                                              std::uint64_t* out) {                 \
        using M = Morton<2, BITS>;                                                  \
        std::array<COORD, 2> a{lo[0], lo[1]}, b{hi[0], hi[1]};                      \
        std::size_t i = 0;                                                          \
        morton::for_each_in_box_zorder<2, BITS>(a, b,                               \
            [&](M m) { out[i++] = m.code(); });                                     \
    }

DEFINE_2D(u32, 32, std::uint32_t)
DEFINE_2D(u16, 16, std::uint16_t)

// ---------------------------------------------------------------------------
// 3D configurations
// ---------------------------------------------------------------------------
#define DEFINE_3D(SUFFIX, BITS, COORD)                                              \
    API void mortonarith_encode3_##SUFFIX(const COORD* x, const COORD* y,           \
                                          const COORD* z, std::uint64_t* out,       \
                                          std::size_t n) {                          \
        enc3<BITS>(x, y, z, out, n);                                                \
    }                                                                               \
    API void mortonarith_decode3_##SUFFIX(const std::uint64_t* code, COORD* x,      \
                                          COORD* y, COORD* z, std::size_t n) {      \
        dec3<BITS>(code, x, y, z, n);                                               \
    }                                                                               \
    API void mortonarith_add3_##SUFFIX(const std::uint64_t* code, std::uint64_t* out, \
                                       std::size_t n, unsigned axis,                 \
                                       std::int64_t k) {                             \
        add3<BITS>(code, out, n, axis, k);                                          \
    }                                                                               \
    API std::uint64_t mortonarith_box_count3_##SUFFIX(const COORD* lo,              \
                                                      const COORD* hi) {            \
        std::uint64_t c = 1;                                                        \
        for (int d = 0; d < 3; ++d) {                                               \
            if (hi[d] < lo[d]) return 0;                                            \
            c *= std::uint64_t(hi[d] - lo[d]) + 1;                                  \
        }                                                                           \
        return c;                                                                   \
    }                                                                               \
    API void mortonarith_box_zorder3_##SUFFIX(const COORD* lo, const COORD* hi,     \
                                              std::uint64_t* out) {                 \
        using M = Morton<3, BITS>;                                                  \
        std::array<COORD, 3> a{lo[0], lo[1], lo[2]}, b{hi[0], hi[1], hi[2]};        \
        std::size_t i = 0;                                                          \
        morton::for_each_in_box_zorder<3, BITS>(a, b,                               \
            [&](M m) { out[i++] = m.code(); });                                     \
    }

DEFINE_3D(u32, 21, std::uint32_t)
DEFINE_3D(u16, 16, std::uint16_t)
