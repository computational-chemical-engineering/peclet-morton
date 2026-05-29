// C ABI shim exposing the morton-arithmetic library to other languages
// (used by the ctypes-based Python package). All functions operate on whole
// arrays so the per-element Python overhead is amortised away and NumPy gets
// near-native throughput.
//
// SPDX-License-Identifier: MIT

#include <cstddef>
#include <cstdint>

#include "morton/iterate.hpp"
#include "morton/morton.hpp"

using morton::Morton;

#define API extern "C" __attribute__((visibility("default")))

// ---------------------------------------------------------------------------
// 2D configurations
// ---------------------------------------------------------------------------
#define DEFINE_2D(SUFFIX, BITS, COORD)                                              \
    API void mortonarith_encode2_##SUFFIX(const COORD* x, const COORD* y,           \
                                          std::uint64_t* out, std::size_t n) {      \
        using M = Morton<2, BITS>;                                                  \
        for (std::size_t i = 0; i < n; ++i)                                         \
            out[i] = M::encode(x[i], y[i]).code();                                  \
    }                                                                               \
    API void mortonarith_decode2_##SUFFIX(const std::uint64_t* code, COORD* x,      \
                                          COORD* y, std::size_t n) {                \
        using M = Morton<2, BITS>;                                                  \
        for (std::size_t i = 0; i < n; ++i) {                                       \
            M m = M::from_code(code[i]);                                            \
            x[i] = COORD(m.get(0));                                                 \
            y[i] = COORD(m.get(1));                                                 \
        }                                                                           \
    }                                                                               \
    API void mortonarith_add2_##SUFFIX(const std::uint64_t* code, std::uint64_t* out, \
                                       std::size_t n, unsigned axis,                 \
                                       std::int64_t k) {                             \
        using M = Morton<2, BITS>;                                                  \
        using C = typename M::coord_type;                                           \
        for (std::size_t i = 0; i < n; ++i) {                                       \
            M m = M::from_code(code[i]);                                            \
            if (k >= 0)                                                             \
                m.add(axis, C(std::uint64_t(k)));                                   \
            else                                                                    \
                m.sub(axis, C(std::uint64_t(-k)));                                  \
            out[i] = m.code();                                                      \
        }                                                                           \
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
        using M = Morton<3, BITS>;                                                  \
        for (std::size_t i = 0; i < n; ++i)                                         \
            out[i] = M::encode(x[i], y[i], z[i]).code();                            \
    }                                                                               \
    API void mortonarith_decode3_##SUFFIX(const std::uint64_t* code, COORD* x,      \
                                          COORD* y, COORD* z, std::size_t n) {      \
        using M = Morton<3, BITS>;                                                  \
        for (std::size_t i = 0; i < n; ++i) {                                       \
            M m = M::from_code(code[i]);                                            \
            x[i] = COORD(m.get(0));                                                 \
            y[i] = COORD(m.get(1));                                                 \
            z[i] = COORD(m.get(2));                                                 \
        }                                                                           \
    }                                                                               \
    API void mortonarith_add3_##SUFFIX(const std::uint64_t* code, std::uint64_t* out, \
                                       std::size_t n, unsigned axis,                 \
                                       std::int64_t k) {                             \
        using M = Morton<3, BITS>;                                                  \
        using C = typename M::coord_type;                                           \
        for (std::size_t i = 0; i < n; ++i) {                                       \
            M m = M::from_code(code[i]);                                            \
            if (k >= 0)                                                             \
                m.add(axis, C(std::uint64_t(k)));                                   \
            else                                                                    \
                m.sub(axis, C(std::uint64_t(-k)));                                  \
            out[i] = m.code();                                                      \
        }                                                                           \
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
