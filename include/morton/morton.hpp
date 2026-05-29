// morton/morton.hpp
//
// A small, header-only C++17 library for Morton (Z-order) codes whose
// distinguishing feature is *arithmetic directly in Morton space*: you can
// increment, decrement, add to or step along a single axis, or find a
// neighbour, in a handful of branchless instructions -- without decoding to
// coordinates and re-encoding.
//
// This makes "walk to the next cell" operations (grid traversal, stencil /
// neighbour access, region fills) essentially free compared with the usual
// decode -> ++ -> encode round trip.
//
// Encoding/decoding uses the BMI2 PDEP/PEXT instructions when available
// (-mbmi2) and falls back to a portable software implementation otherwise.
//
// Constraints: Dim * Bits <= 64. (1, 2 or 3 dimensions are the common cases,
// but any Dim that fits is supported.)
//
// SPDX-License-Identifier: MIT

#ifndef MORTON_MORTON_HPP
#define MORTON_MORTON_HPP

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

#if defined(__BMI2__)
#include <immintrin.h>
#endif

namespace morton {

namespace detail {

// Smallest unsigned integer type that holds at least NBits bits.
template <unsigned NBits>
struct uint_for {
    using type = std::conditional_t<
        (NBits <= 8), std::uint8_t,
        std::conditional_t<
            (NBits <= 16), std::uint16_t,
            std::conditional_t<(NBits <= 32), std::uint32_t, std::uint64_t>>>;
};
template <unsigned NBits>
using uint_for_t = typename uint_for<NBits>::type;

// Build the mask selecting the bits that belong to axis `d` in a code with
// `Dim` interleaved axes of `Bits` bits each: positions d, d+Dim, d+2*Dim, ...
template <unsigned Dim, unsigned Bits, typename Code>
constexpr Code axis_mask(unsigned d) {
    Code m = 0;
    for (unsigned i = 0; i < Bits; ++i)
        m |= (Code(1) << (i * Dim + d));
    return m;
}

// Software deposit: spread the low `Bits` bits of x across positions
// d, d+Dim, ... (the portable fallback for PDEP). Unrolls to a chain of
// shifts/ors because Bits is a compile-time constant.
template <unsigned Dim, unsigned Bits, typename Code, typename Coord>
constexpr Code spread_sw(Coord x) {
    Code r = 0;
    for (unsigned i = 0; i < Bits; ++i)
        r |= Code(Code(x >> i) & Code(1)) << (i * Dim);
    return r;
}

// Software extract: inverse of spread_sw (the portable fallback for PEXT).
template <unsigned Dim, unsigned Bits, typename Code, typename Coord>
constexpr Coord compact_sw(Code c) {
    Coord r = 0;
    for (unsigned i = 0; i < Bits; ++i)
        r |= Coord(Coord(c >> (i * Dim)) & Coord(1)) << i;
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
    static_assert(Dim * Bits <= 64, "Dim * Bits must be <= 64");

public:
    static constexpr unsigned dimensions = Dim;
    static constexpr unsigned bits_per_axis = Bits;
    static constexpr unsigned code_bits = Dim * Bits;

    using code_type = detail::uint_for_t<Dim * Bits>;
    using coord_type = detail::uint_for_t<Bits>;

    /// Mask of all valid code bits (handles code_bits == 64 without UB).
    static constexpr code_type field_mask =
        (code_bits >= std::numeric_limits<code_type>::digits)
            ? code_type(~code_type(0))
            : code_type((code_type(1) << code_bits) - 1);

    /// Largest representable coordinate value per axis.
    static constexpr coord_type coord_max =
        (Bits >= std::numeric_limits<coord_type>::digits)
            ? coord_type(~coord_type(0))
            : coord_type((coord_type(1) << Bits) - 1);

    // ---- construction -----------------------------------------------------

    constexpr Morton() noexcept : code_(0) {}

    /// Wrap a raw, already-interleaved code value.
    static constexpr Morton from_code(code_type raw) noexcept {
        Morton m;
        m.code_ = raw & field_mask;
        return m;
    }

    /// Encode `Dim` coordinates into a Morton code.
    template <typename... Cs,
              typename = std::enable_if_t<sizeof...(Cs) == Dim>>
    static Morton encode(Cs... coords) {
        coord_type tmp[Dim] = {static_cast<coord_type>(coords)...};
        Morton m;
        for (unsigned d = 0; d < Dim; ++d)
            m.code_ |= deposit(tmp[d], d);
        return m;
    }

    /// Encode from an array of coordinates.
    static Morton encode(const std::array<coord_type, Dim>& c) {
        Morton m;
        for (unsigned d = 0; d < Dim; ++d)
            m.code_ |= deposit(c[d], d);
        return m;
    }

    // ---- access -----------------------------------------------------------

    constexpr code_type code() const noexcept { return code_; }

    /// Decode coordinate of a single axis.
    coord_type get(unsigned d) const { return extract(code_, d); }

    /// Decode all coordinates.
    std::array<coord_type, Dim> decode() const {
        std::array<coord_type, Dim> out{};
        for (unsigned d = 0; d < Dim; ++d)
            out[d] = extract(code_, d);
        return out;
    }

    /// Replace one axis' coordinate (other axes untouched). Wraps mod 2^Bits.
    void set(unsigned d, coord_type value) {
        code_ = (code_ & keep_mask(d)) | deposit(value, d);
    }

    // ---- the headline: arithmetic in Morton space -------------------------
    //
    // Each of these touches only the bits of one axis and leaves the other
    // axes' interleaved bits exactly where they are -- no decode/encode.

    /// Add `k` to axis `d` (wraps mod 2^Bits). O(1), branchless.
    void add(unsigned d, coord_type k) {
        const code_type M = axis_mask(d);
        // Fill the non-axis bits with 1s so carries ripple across the gaps,
        // add the dilated increment, then keep only this axis' result.
        code_type s = (code_ | ~M) + deposit(k, d);
        code_ = (s & M) | (code_ & M_complement(d));
    }

    /// Subtract `k` from axis `d` (wraps mod 2^Bits). O(1), branchless.
    void sub(unsigned d, coord_type k) {
        const code_type M = axis_mask(d);
        // Non-axis bits are 0 here, so borrows ripple across the gaps.
        code_type s = (code_ & M) - deposit(k, d);
        code_ = (s & M) | (code_ & M_complement(d));
    }

    /// Increment axis `d` by one. O(1).
    void inc(unsigned d) {
        const code_type M = axis_mask(d);
        code_type s = (code_ | ~M) + lsb(d);
        code_ = (s & M) | (code_ & M_complement(d));
    }

    /// Decrement axis `d` by one. O(1).
    void dec(unsigned d) {
        const code_type M = axis_mask(d);
        code_type s = (code_ & M) - lsb(d);
        code_ = (s & M) | (code_ & M_complement(d));
    }

    /// Morton code of the neighbour one step along `±` axis `d` (wraps).
    Morton neighbor(unsigned d, int dir) const {
        Morton m = *this;
        if (dir >= 0)
            m.inc(d);
        else
            m.dec(d);
        return m;
    }

    // ---- Z-order successor / predecessor ----------------------------------
    //
    // The Morton code *is* the Z-order index, so the next/previous cell in
    // Z-order is just ++ / -- on the integer.

    Morton& operator++() noexcept {
        code_ = (code_ + 1) & field_mask;
        return *this;
    }
    Morton operator++(int) noexcept {
        Morton t = *this;
        ++*this;
        return t;
    }
    Morton& operator--() noexcept {
        code_ = (code_ - 1) & field_mask;
        return *this;
    }
    Morton operator--(int) noexcept {
        Morton t = *this;
        --*this;
        return t;
    }

    // ---- comparison (Z-order) ---------------------------------------------

    friend constexpr bool operator==(Morton a, Morton b) noexcept {
        return a.code_ == b.code_;
    }
    friend constexpr bool operator!=(Morton a, Morton b) noexcept {
        return a.code_ != b.code_;
    }
    friend constexpr bool operator<(Morton a, Morton b) noexcept {
        return a.code_ < b.code_;
    }
    friend constexpr bool operator<=(Morton a, Morton b) noexcept {
        return a.code_ <= b.code_;
    }
    friend constexpr bool operator>(Morton a, Morton b) noexcept {
        return a.code_ > b.code_;
    }
    friend constexpr bool operator>=(Morton a, Morton b) noexcept {
        return a.code_ >= b.code_;
    }

    // ---- low-level deposit/extract (public for the bindings/benchmarks) ---

    static code_type deposit(coord_type x, unsigned d) {
#if defined(__BMI2__)
        return code_type(_pdep_u64(std::uint64_t(x), std::uint64_t(axis_mask(d))));
#else
        return detail::spread_sw<Dim, Bits, code_type, coord_type>(x) << d;
#endif
    }

    static coord_type extract(code_type c, unsigned d) {
#if defined(__BMI2__)
        return coord_type(_pext_u64(std::uint64_t(c), std::uint64_t(axis_mask(d))));
#else
        return detail::compact_sw<Dim, Bits, code_type, coord_type>(c >> d);
#endif
    }

    /// Mask of bits belonging to axis `d`.
    static code_type axis_mask(unsigned d) {
        return mask_table()[d];
    }

private:
    static code_type M_complement(unsigned d) {
        return field_mask & ~axis_mask(d);
    }
    static code_type keep_mask(unsigned d) { return M_complement(d); }

    // Lowest bit of axis d.
    static constexpr code_type lsb(unsigned d) { return code_type(1) << d; }

    // Per-axis masks, computed once.
    static const std::array<code_type, Dim>& mask_table() {
        static const std::array<code_type, Dim> t = [] {
            std::array<code_type, Dim> a{};
            for (unsigned d = 0; d < Dim; ++d)
                a[d] = detail::axis_mask<Dim, Bits, code_type>(d);
            return a;
        }();
        return t;
    }

    code_type code_;
};

// Convenient aliases for the common cases.
using Morton2D32 = Morton<2, 32>;  // 2D, 32 bits/axis  -> 64-bit code
using Morton2D16 = Morton<2, 16>;  // 2D, 16 bits/axis  -> 32-bit code
using Morton3D21 = Morton<3, 21>;  // 3D, 21 bits/axis  -> 63-bit code
using Morton3D16 = Morton<3, 16>;  // 3D, 16 bits/axis  -> 48-bit code

}  // namespace morton

#endif  // MORTON_MORTON_HPP
