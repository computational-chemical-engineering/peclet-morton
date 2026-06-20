/// @file wide_uint.hpp
/// @brief Minimal fixed-width unsigned integer backing Morton codes wider than 128 bits.
///
/// `detail::wide_uint<W>` is `W` 64-bit words (little-endian) providing exactly the operators
/// `Morton<Dim,Bits>` uses: `+ - & | ^ ~ << >>` and comparisons. This lets the same Morton
/// implementation work for codes wider than 128 bits, where no built-in integer exists. It is
/// intentionally small and `constexpr`; it is **not** a general bignum and should not be tuned for
/// speed unless it becomes a hot path.
///
/// SPDX-License-Identifier: MIT

#ifndef MORTON_WIDE_UINT_HPP
#define MORTON_WIDE_UINT_HPP

#include <array>
#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace morton {
namespace detail {

/// Fixed-width unsigned integer of @tparam W 64-bit words, little-endian (`w[0]` least significant).
/// Selected automatically by `uint_for` when `Dim*Bits` exceeds the largest built-in integer.
template <std::size_t W>
struct wide_uint {
    static_assert(W >= 2, "use a built-in integer for <= 64 bits");
    std::array<std::uint64_t, W> w{};  ///< Words, little-endian: w[0] is least significant.

    constexpr wide_uint() = default;

    // From any built-in integral value up to 64 bits (low word).
    template <typename T,
              typename = std::enable_if_t<std::is_integral_v<T> && (sizeof(T) <= 8)>>
    constexpr wide_uint(T v) : w{} {
        w[0] = std::uint64_t(static_cast<std::make_unsigned_t<T>>(v));
    }

#if defined(__SIZEOF_INT128__)
    // From a 128-bit built-in (fills the low two words).
    constexpr wide_uint(unsigned __int128 v) : w{} {
        w[0] = std::uint64_t(v);
        if constexpr (W >= 2) w[1] = std::uint64_t(v >> 64);
    }
#endif

    // From a wide_uint of a different width (truncate or zero-extend).
    template <std::size_t W2, typename = std::enable_if_t<W2 != W>>
    constexpr explicit wide_uint(const wide_uint<W2>& o) : w{} {
        for (std::size_t i = 0; i < W && i < W2; ++i) w[i] = o.w[i];
    }

    explicit constexpr operator std::uint64_t() const { return w[0]; }
#if defined(__SIZEOF_INT128__)
    explicit constexpr operator unsigned __int128() const {
        unsigned __int128 r = w[0];
        if constexpr (W >= 2) r |= (unsigned __int128)(w[1]) << 64;
        return r;
    }
#endif
    explicit constexpr operator bool() const {
        for (std::size_t i = 0; i < W; ++i)
            if (w[i]) return true;
        return false;
    }

    // ---- bitwise ----------------------------------------------------------
    constexpr wide_uint operator~() const {
        wide_uint r;
        for (std::size_t i = 0; i < W; ++i) r.w[i] = ~w[i];
        return r;
    }
    constexpr wide_uint& operator&=(const wide_uint& o) {
        for (std::size_t i = 0; i < W; ++i) w[i] &= o.w[i];
        return *this;
    }
    constexpr wide_uint& operator|=(const wide_uint& o) {
        for (std::size_t i = 0; i < W; ++i) w[i] |= o.w[i];
        return *this;
    }
    constexpr wide_uint& operator^=(const wide_uint& o) {
        for (std::size_t i = 0; i < W; ++i) w[i] ^= o.w[i];
        return *this;
    }
    friend constexpr wide_uint operator&(wide_uint a, const wide_uint& b) { return a &= b; }
    friend constexpr wide_uint operator|(wide_uint a, const wide_uint& b) { return a |= b; }
    friend constexpr wide_uint operator^(wide_uint a, const wide_uint& b) { return a ^= b; }

    // ---- shifts -----------------------------------------------------------
    constexpr wide_uint operator<<(unsigned s) const {
        wide_uint r;
        const unsigned ws = s / 64, bs = s % 64;
        for (int i = int(W) - 1; i >= 0; --i) {
            const int src = i - int(ws);
            std::uint64_t v = 0;
            if (src >= 0) {
                v = w[std::size_t(src)] << bs;
                if (bs && src - 1 >= 0) v |= w[std::size_t(src - 1)] >> (64 - bs);
            }
            r.w[std::size_t(i)] = v;
        }
        return r;
    }
    constexpr wide_uint operator>>(unsigned s) const {
        wide_uint r;
        const unsigned ws = s / 64, bs = s % 64;
        for (std::size_t i = 0; i < W; ++i) {
            const std::size_t src = i + ws;
            std::uint64_t v = 0;
            if (src < W) {
                v = w[src] >> bs;
                if (bs && src + 1 < W) v |= w[src + 1] << (64 - bs);
            }
            r.w[i] = v;
        }
        return r;
    }
    constexpr wide_uint& operator<<=(unsigned s) { return *this = (*this << s); }
    constexpr wide_uint& operator>>=(unsigned s) { return *this = (*this >> s); }

    // ---- add / sub (ripple) -----------------------------------------------
    friend constexpr wide_uint operator+(const wide_uint& a, const wide_uint& b) {
        wide_uint r;
        std::uint64_t carry = 0;
        for (std::size_t i = 0; i < W; ++i) {
            std::uint64_t s1 = a.w[i] + b.w[i];
            std::uint64_t c1 = (s1 < a.w[i]) ? 1u : 0u;
            std::uint64_t s2 = s1 + carry;
            std::uint64_t c2 = (s2 < s1) ? 1u : 0u;
            r.w[i] = s2;
            carry = c1 | c2;
        }
        return r;
    }
    friend constexpr wide_uint operator-(const wide_uint& a, const wide_uint& b) {
        wide_uint r;
        std::uint64_t borrow = 0;
        for (std::size_t i = 0; i < W; ++i) {
            std::uint64_t d1 = a.w[i] - b.w[i];
            std::uint64_t b1 = (a.w[i] < b.w[i]) ? 1u : 0u;
            std::uint64_t d2 = d1 - borrow;
            std::uint64_t b2 = (d1 < borrow) ? 1u : 0u;
            r.w[i] = d2;
            borrow = b1 | b2;
        }
        return r;
    }
    constexpr wide_uint& operator++() { return *this = (*this + wide_uint(1)); }
    constexpr wide_uint& operator--() { return *this = (*this - wide_uint(1)); }

    // ---- comparison -------------------------------------------------------
    friend constexpr bool operator==(const wide_uint& a, const wide_uint& b) {
        for (std::size_t i = 0; i < W; ++i)
            if (a.w[i] != b.w[i]) return false;
        return true;
    }
    friend constexpr bool operator!=(const wide_uint& a, const wide_uint& b) { return !(a == b); }
    friend constexpr bool operator<(const wide_uint& a, const wide_uint& b) {
        for (std::size_t i = W; i-- > 0;) {
            if (a.w[i] != b.w[i]) return a.w[i] < b.w[i];
        }
        return false;
    }
    friend constexpr bool operator>(const wide_uint& a, const wide_uint& b) { return b < a; }
    friend constexpr bool operator<=(const wide_uint& a, const wide_uint& b) { return !(b < a); }
    friend constexpr bool operator>=(const wide_uint& a, const wide_uint& b) { return !(a < b); }
};

}  // namespace detail
}  // namespace morton

#endif  // MORTON_WIDE_UINT_HPP
