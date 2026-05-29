// morton/morton.hpp
//
// A small, header-only C++17 library for Morton (Z-order) codes whose
// distinguishing feature is *arithmetic directly in Morton space*: you can
// increment, decrement, add to or step along a single axis, or find a
// neighbour, in a handful of branchless instructions -- without decoding to
// coordinates and re-encoding.
//
// Encoding/decoding uses the BMI2 PDEP/PEXT instructions when available
// (-mbmi2, codes up to 64 bits) and a portable software implementation
// otherwise. The software path is constexpr, so encode/decode/arithmetic can
// run at compile time (e.g. to build lookup tables).
//
// Codes up to 64 bits are stored in a built-in unsigned integer; codes from 65
// to 128 bits use __uint128_t where the compiler provides it (GCC/Clang), so
// 3D 32-bit (96 bits) and 2D 64-bit (128 bits) work too.
//
// SPDX-License-Identifier: MIT

#ifndef MORTON_MORTON_HPP
#define MORTON_MORTON_HPP

#include <array>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <type_traits>

#include "morton/wide_uint.hpp"

#if defined(__BMI2__)
#include <immintrin.h>
#endif

#if defined(__SIZEOF_INT128__)
#define MORTON_HAS_INT128 1
#endif

// Maximum supported code width (Dim * Bits). Codes <= 64 use a built-in
// integer, <= 128 use __uint128_t (where available), and wider codes use the
// software wide_uint<W> backend.
#ifndef MORTON_MAX_BITS
#define MORTON_MAX_BITS 256
#endif

// Mark functions callable from CUDA/HIP device code. Expands to nothing for an
// ordinary host C++ build, so the CPU library is completely unchanged.
#if defined(__CUDACC__) || defined(__HIPCC__)
#define MORTON_HD __host__ __device__
#else
#define MORTON_HD
#endif

namespace morton {

namespace detail {

#if defined(MORTON_HAS_INT128)
using uint128_t = unsigned __int128;
#endif

// Smallest unsigned type that holds at least NBits bits: a built-in where one
// exists, otherwise a wide_uint of the right number of 64-bit words.
template <unsigned NBits>
struct uint_for {
#if defined(MORTON_HAS_INT128)
    static constexpr unsigned builtin_max = 128;
#else
    static constexpr unsigned builtin_max = 64;
#endif
    using builtin = std::conditional_t<
        (NBits <= 8), std::uint8_t,
        std::conditional_t<
            (NBits <= 16), std::uint16_t,
            std::conditional_t<(NBits <= 32), std::uint32_t,
#if defined(MORTON_HAS_INT128)
                               std::conditional_t<(NBits <= 64), std::uint64_t, uint128_t>
#else
                               std::uint64_t
#endif
                               >>>;
    using type = std::conditional_t<(NBits <= builtin_max), builtin,
                                     wide_uint<(NBits + 63) / 64>>;
};
template <unsigned NBits>
using uint_for_t = typename uint_for<NBits>::type;

// True during constant evaluation. Lets a constexpr function avoid the
// (non-constexpr) BMI2 intrinsics when evaluated at compile time.
MORTON_HD constexpr bool is_consteval() noexcept {
#if defined(__cpp_if_consteval)
    if consteval {
        return true;
    } else {
        return false;
    }
#elif defined(__GNUC__) || defined(__clang__)
    return __builtin_is_constant_evaluated();
#else
    return false;
#endif
}

// 3^n, as a compile-time value (size of a Moore neighbourhood including self).
MORTON_HD constexpr unsigned pow3(unsigned n) {
    unsigned r = 1;
    for (unsigned i = 0; i < n; ++i) r *= 3;
    return r;
}

// Build the mask selecting the bits that belong to axis `d` in a code with
// `Dim` interleaved axes of `Bits` bits each: positions d, d+Dim, d+2*Dim, ...
template <unsigned Dim, unsigned Bits, typename Code>
MORTON_HD constexpr Code axis_mask(unsigned d) {
    Code m = 0;
    for (unsigned i = 0; i < Bits; ++i)
        m |= (Code(1) << (i * Dim + d));
    return m;
}

// Software deposit: spread the low `Bits` bits of x across positions
// d, d+Dim, ... (the portable fallback for PDEP). Works for any code width.
template <unsigned Dim, unsigned Bits, typename Code, typename Coord>
MORTON_HD constexpr Code spread_sw(Coord x, unsigned d) {
    Code r = 0;
    for (unsigned i = 0; i < Bits; ++i)
        r |= Code(Code(x >> i) & Code(1)) << (i * Dim + d);
    return r;
}

// Software extract: inverse of spread_sw (the portable fallback for PEXT).
template <unsigned Dim, unsigned Bits, typename Code, typename Coord>
MORTON_HD constexpr Coord compact_sw(Code c, unsigned d) {
    Coord r = 0;
    for (unsigned i = 0; i < Bits; ++i)
        r |= Coord(Coord(c >> (i * Dim + d)) & Coord(1)) << i;
    return r;
}

}  // namespace detail

/// A Morton (Z-order) code interleaving `Dim` unsigned coordinates of `Bits`
/// bits each. The whole code occupies `Dim * Bits` bits and is stored in the
/// smallest unsigned integer that fits.
template <unsigned Dim, unsigned Bits>
class Morton {
    static_assert(Dim >= 1, "Dim must be >= 1");
    static_assert(Bits >= 1, "Bits must be >= 1");
    static_assert(Dim * Bits <= MORTON_MAX_BITS,
                  "Dim * Bits exceeds MORTON_MAX_BITS (raise it if you need wider codes)");

public:
    static constexpr unsigned dimensions = Dim;
    static constexpr unsigned bits_per_axis = Bits;
    static constexpr unsigned code_bits = Dim * Bits;

    using code_type = detail::uint_for_t<Dim * Bits>;
    using coord_type = detail::uint_for_t<Bits>;

private:
    static constexpr unsigned code_type_bits = sizeof(code_type) * CHAR_BIT;
    static constexpr unsigned coord_type_bits = sizeof(coord_type) * CHAR_BIT;

public:
    /// Mask of all valid code bits (handles code_bits == width without UB).
    static constexpr code_type field_mask =
        (code_bits >= code_type_bits) ? code_type(~code_type(0))
                                      : code_type((code_type(1) << code_bits) - 1);

    /// Largest representable coordinate value per axis.
    static constexpr coord_type coord_max =
        (Bits >= coord_type_bits) ? coord_type(~coord_type(0))
                                  : coord_type((coord_type(1) << Bits) - 1);

    static constexpr unsigned octants = unsigned(1) << Dim;  // children per node

    // ---- construction -----------------------------------------------------

    MORTON_HD constexpr Morton() noexcept : code_(0) {}

    /// Wrap a raw, already-interleaved code value.
    MORTON_HD static constexpr Morton from_code(code_type raw) noexcept {
        Morton m;
        m.code_ = raw & field_mask;
        return m;
    }

    /// Encode `Dim` coordinates into a Morton code.
    template <typename... Cs, typename = std::enable_if_t<sizeof...(Cs) == Dim>>
    MORTON_HD static constexpr Morton encode(Cs... coords) {
        coord_type tmp[Dim] = {static_cast<coord_type>(coords)...};
        Morton m;
        for (unsigned d = 0; d < Dim; ++d)
            m.code_ |= deposit(tmp[d], d);
        return m;
    }

    /// Encode from an array of coordinates.
    MORTON_HD static constexpr Morton encode(const std::array<coord_type, Dim>& c) {
        Morton m;
        for (unsigned d = 0; d < Dim; ++d)
            m.code_ |= deposit(c[d], d);
        return m;
    }

    // ---- access -----------------------------------------------------------

    MORTON_HD constexpr code_type code() const noexcept { return code_; }

    /// Decode coordinate of a single axis.
    MORTON_HD constexpr coord_type get(unsigned d) const { return extract(code_, d); }

    /// Decode all coordinates.
    MORTON_HD constexpr std::array<coord_type, Dim> decode() const {
        std::array<coord_type, Dim> out{};
        for (unsigned d = 0; d < Dim; ++d)
            out[d] = extract(code_, d);
        return out;
    }

    /// Replace one axis' coordinate (other axes untouched). Wraps mod 2^Bits.
    MORTON_HD constexpr void set(unsigned d, coord_type value) {
        code_ = (code_ & keep_mask(d)) | deposit(value, d);
    }

    // ---- the headline: arithmetic in Morton space -------------------------
    //
    // Each of these touches only the bits of one axis and leaves the other
    // axes' interleaved bits exactly where they are -- no decode/encode.
    // All wrap modulo 2^Bits on the affected axis.

    /// Add `k` to axis `d` (wraps). O(1), branchless.
    MORTON_HD constexpr void add(unsigned d, coord_type k) {
        const code_type M = axis_mask(d);
        code_type s = (code_ | ~M) + deposit(k, d);
        code_ = (s & M) | (code_ & M_complement(d));
    }

    /// Subtract `k` from axis `d` (wraps). O(1), branchless.
    MORTON_HD constexpr void sub(unsigned d, coord_type k) {
        const code_type M = axis_mask(d);
        code_type s = (code_ & M) - deposit(k, d);
        code_ = (s & M) | (code_ & M_complement(d));
    }

    /// Increment axis `d` by one. O(1).
    MORTON_HD constexpr void inc(unsigned d) {
        const code_type M = axis_mask(d);
        code_type s = (code_ | ~M) + lsb(d);
        code_ = (s & M) | (code_ & M_complement(d));
    }

    /// Decrement axis `d` by one. O(1).
    MORTON_HD constexpr void dec(unsigned d) {
        const code_type M = axis_mask(d);
        code_type s = (code_ & M) - lsb(d);
        code_ = (s & M) | (code_ & M_complement(d));
    }

    /// Morton code of the neighbour one step along `±` axis `d` (wraps).
    MORTON_HD constexpr Morton neighbor(unsigned d, int dir) const {
        Morton m = *this;
        if (dir >= 0)
            m.inc(d);
        else
            m.dec(d);
        return m;
    }

    // ---- saturating / checked arithmetic (do not wrap) --------------------

    /// Add `k` to axis `d`, clamping at coord_max instead of wrapping.
    MORTON_HD constexpr void add_sat(unsigned d, coord_type k) {
        coord_type cur = get(d);
        coord_type room = coord_type(coord_max - cur);
        set(d, (k > room) ? coord_max : coord_type(cur + k));
    }

    /// Subtract `k` from axis `d`, clamping at 0 instead of wrapping.
    MORTON_HD constexpr void sub_sat(unsigned d, coord_type k) {
        coord_type cur = get(d);
        set(d, (k > cur) ? coord_type(0) : coord_type(cur - k));
    }

    /// Add `k` to axis `d` only if it does not overflow past coord_max.
    /// Returns false (and leaves the code unchanged) if it would.
    MORTON_HD constexpr bool try_add(unsigned d, coord_type k) {
        coord_type cur = get(d);
        if (k > coord_type(coord_max - cur)) return false;
        set(d, coord_type(cur + k));
        return true;
    }

    /// Subtract `k` from axis `d` only if it does not underflow past 0.
    MORTON_HD constexpr bool try_sub(unsigned d, coord_type k) {
        coord_type cur = get(d);
        if (k > cur) return false;
        set(d, coord_type(cur - k));
        return true;
    }

    // ---- neighbour sets ----------------------------------------------------

    /// The 2*Dim von Neumann (face) neighbours, ordered
    /// {axis0-, axis0+, axis1-, axis1+, ...}. Wraps at the grid edge.
    MORTON_HD constexpr std::array<Morton, 2 * Dim> face_neighbors() const {
        std::array<Morton, 2 * Dim> out{};
        for (unsigned d = 0; d < Dim; ++d) {
            out[2 * d] = neighbor(d, -1);
            out[2 * d + 1] = neighbor(d, +1);
        }
        return out;
    }

    /// The 3^Dim - 1 Moore neighbours (all cells differing by -1/0/+1 on each
    /// axis, excluding self), in odometer order. Wraps at the grid edge.
    MORTON_HD constexpr std::array<Morton, detail::pow3(Dim) - 1> all_neighbors() const {
        std::array<Morton, detail::pow3(Dim) - 1> out{};
        unsigned n = 0;
        for (unsigned i = 0; i < detail::pow3(Dim); ++i) {
            int off[Dim];
            unsigned t = i;
            bool zero = true;
            for (unsigned d = 0; d < Dim; ++d) {
                off[d] = int(t % 3) - 1;
                if (off[d] != 0) zero = false;
                t /= 3;
            }
            if (zero) continue;  // skip self
            Morton m = *this;
            for (unsigned d = 0; d < Dim; ++d) {
                if (off[d] > 0)
                    m.inc(d);
                else if (off[d] < 0)
                    m.dec(d);
            }
            out[n++] = m;
        }
        return out;
    }

    // ---- octree / hierarchy navigation ------------------------------------
    //
    // A cell at `level` covers a 2^level block per axis and has its low
    // level*Dim code bits zero. These helpers express that hierarchy directly.

    /// Ancestor cell origin at the given `level` (clears the low level*Dim
    /// bits). level 0 is `*this`.
    MORTON_HD constexpr Morton ancestor(unsigned level) const {
        if (level * Dim >= code_bits) return Morton{};
        code_type lowmask = (code_type(1) << (level * Dim)) - 1;
        return from_code(code_ & ~lowmask);
    }

    /// Index (0 .. 2^Dim-1) of this cell within its parent at `level+1`,
    /// i.e. the Dim interleaved bits just above the level boundary.
    MORTON_HD constexpr unsigned child_index(unsigned level) const {
        return unsigned((code_ >> (level * Dim)) & (octants - 1));
    }

    /// The `octant`-th child origin of a cell that is an ancestor at `level`
    /// (i.e. set the Dim bits at position (level-1)*Dim). Requires level >= 1.
    MORTON_HD constexpr Morton child(unsigned level, unsigned octant) const {
        unsigned shift = (level - 1) * Dim;
        code_type cleared = code_ & ~(code_type(octants - 1) << shift);
        return from_code(cleared | (code_type(octant) << shift));
    }

    // ---- Z-order successor / predecessor ----------------------------------

    MORTON_HD constexpr Morton& operator++() noexcept {
        code_ = (code_ + 1) & field_mask;
        return *this;
    }
    MORTON_HD constexpr Morton operator++(int) noexcept {
        Morton t = *this;
        ++*this;
        return t;
    }
    MORTON_HD constexpr Morton& operator--() noexcept {
        code_ = (code_ - 1) & field_mask;
        return *this;
    }
    MORTON_HD constexpr Morton operator--(int) noexcept {
        Morton t = *this;
        --*this;
        return t;
    }

    // ---- comparison (Z-order) ---------------------------------------------

    friend constexpr bool operator==(Morton a, Morton b) noexcept { return a.code_ == b.code_; }
    friend constexpr bool operator!=(Morton a, Morton b) noexcept { return a.code_ != b.code_; }
    friend MORTON_HD constexpr bool operator<(Morton a, Morton b) noexcept { return a.code_ < b.code_; }
    friend constexpr bool operator<=(Morton a, Morton b) noexcept { return a.code_ <= b.code_; }
    friend MORTON_HD constexpr bool operator>(Morton a, Morton b) noexcept { return a.code_ > b.code_; }
    friend constexpr bool operator>=(Morton a, Morton b) noexcept { return a.code_ >= b.code_; }

    // ---- low-level deposit/extract ----------------------------------------

    MORTON_HD static constexpr code_type deposit(coord_type x, unsigned d) {
        if constexpr (code_bits <= 64) {
            // PDEP is a host x86 instruction; never emit it in CUDA device code
            // (where __CUDA_ARCH__ is defined), even if the host was -mbmi2.
#if defined(__BMI2__) && !defined(__CUDA_ARCH__)
            if (!detail::is_consteval())
                return code_type(_pdep_u64(std::uint64_t(x), std::uint64_t(axis_mask(d))));
#endif
        }
        return detail::spread_sw<Dim, Bits, code_type, coord_type>(x, d);
    }

    MORTON_HD static constexpr coord_type extract(code_type c, unsigned d) {
        if constexpr (code_bits <= 64) {
#if defined(__BMI2__) && !defined(__CUDA_ARCH__)
            if (!detail::is_consteval())
                return coord_type(_pext_u64(std::uint64_t(c), std::uint64_t(axis_mask(d))));
#endif
        }
        return detail::compact_sw<Dim, Bits, code_type, coord_type>(c, d);
    }

    /// Mask of bits belonging to axis `d`.
    MORTON_HD static constexpr code_type axis_mask(unsigned d) {
#if defined(__CUDA_ARCH__)
        // On device, compute the mask (avoids referencing the host static
        // array from device code); it is a short compile-time-sized loop.
        return detail::axis_mask<Dim, Bits, code_type>(d);
#else
        return axis_masks_[d];
#endif
    }

private:
    MORTON_HD static constexpr code_type M_complement(unsigned d) { return field_mask & ~axis_mask(d); }
    MORTON_HD static constexpr code_type keep_mask(unsigned d) { return M_complement(d); }
    MORTON_HD static constexpr code_type lsb(unsigned d) { return code_type(1) << d; }

    MORTON_HD static constexpr std::array<code_type, Dim> make_masks() {
        std::array<code_type, Dim> a{};
        for (unsigned d = 0; d < Dim; ++d)
            a[d] = detail::axis_mask<Dim, Bits, code_type>(d);
        return a;
    }
    static constexpr std::array<code_type, Dim> axis_masks_ = make_masks();

    code_type code_;
};

// Convenient aliases for the common cases.
using Morton2D32 = Morton<2, 32>;  // 2D, 32 bits/axis  -> 64-bit code
using Morton2D16 = Morton<2, 16>;  // 2D, 16 bits/axis  -> 32-bit code
using Morton3D21 = Morton<3, 21>;  // 3D, 21 bits/axis  -> 63-bit code
using Morton3D16 = Morton<3, 16>;  // 3D, 16 bits/axis  -> 48-bit code
#if defined(MORTON_HAS_INT128)
using Morton3D32 = Morton<3, 32>;  // 3D, 32 bits/axis  -> 96-bit code
using Morton2D64 = Morton<2, 64>;  // 2D, 64 bits/axis  -> 128-bit code
#endif

}  // namespace morton

#endif  // MORTON_MORTON_HPP
