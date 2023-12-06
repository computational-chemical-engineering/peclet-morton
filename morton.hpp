#ifndef MORTON_HPP
#define MORTON_HPP

#include <array>
#include <limits>
#include <cstring>
#include <utility>  
#include "bitarray.hpp"

using std::array;

template <size_t Dim, size_t N>
class Morton : public BitArrayBase<Morton<Dim, N>, Dim * N>
{
public:
    // Bring constructors from BitArrayBase into BitArray
    using BitArrayBase<Morton<Dim, N>, Dim * N>::BitArrayBase;

    // Constructor to convert from BitArrayBase to Morton
    Morton(const BitArrayBase<Morton<Dim, N>, Dim * N> &other)
        : BitArrayBase<Morton<Dim, N>, Dim * N>(other) {}

    Morton<Dim, N> &operator<<=(size_t pos)
    {
        BitArrayBase<Morton<Dim, N>, Dim * N>::operator<<=(Dim * pos);
        return *this;
    }

    Morton<Dim, N> &operator>>=(size_t pos)
    {
        BitArrayBase<Morton<Dim, N>, Dim * N>::operator>>=(Dim * pos);
        return *this;
    }

    Morton<Dim, N> operator+=(const Morton<Dim, N> &other)
    {
        Morton<Dim, N> carry = *this & other; // Initial carry
        *this ^= other;                       // Initial sum (without carry)

        // Propagate carry without creating additional Morton objects
        while (carry.any())
        {
            Morton<Dim, N> shiftedCarry = carry;
            shiftedCarry <<= 1;
//            Morton<Dim, N> newCarry = *this & shiftedCarry;
            carry = *this & shiftedCarry;
            *this ^= shiftedCarry;
//            carry = newCarry;
        }

        return *this;
    }

    Morton &operator-=(const Morton &other)
    {
        Morton<Dim, N> ones = static_cast<Morton<Dim, N>>((static_cast<size_t>(1) << Dim) - 1);
        Morton<Dim, N> complement = (~other + ones);
        operator+=(complement);
        return *this;
    }

    bool test(size_t dir, size_t pos) const
    {
        return BitArrayBase<Morton<Dim, N>, Dim * N>::test(Dim * pos + dir);
    }

    Morton &set(size_t dir, size_t pos, bool value = true)
    {
        BitArrayBase<Morton<Dim, N>, Dim * N>::set(Dim * pos + dir, value);
        return *this;
    }

    Morton &reset(size_t dir, size_t pos)
    {
        BitArrayBase<Morton<Dim, N>, Dim * N>::reset(Dim * pos + dir);
        return *this;
    }

    size_t countr_zero() const
    {
        return (BitArrayBase<Morton<Dim, N>, Dim * N>::countr_zero()/Dim);
    }

    static constexpr size_t szInBytes = (Dim * N + 7) / 8;
    static constexpr size_t szBitsetInBytes = (N + 7) / 8;
};

template <size_t Dim, size_t N>
Morton<Dim, N> operator<<(const Morton<Dim, N> &m, size_t shift)
{
    Morton<Dim, N> result = m; // Create a copy of the input Morton object
    result <<= shift;          // Apply the left shift
    return result;             // Return the modified copy
}

template <size_t Dim, size_t N>
Morton<Dim, N> operator>>(const Morton<Dim, N> &m, size_t shift)
{
    Morton<Dim, N> result = m; // Create a copy of the input Morton object
    result >>= shift;          // Apply the left shift
    return result;             // Return the modified copy
}

template <size_t Dim, size_t N>
Morton<Dim, N> operator+(const Morton<Dim, N> &lhs, const Morton<Dim, N> &rhs)
{
    Morton<Dim, N> carry = lhs & rhs;  // Initial carry
    Morton<Dim, N> result = lhs ^ rhs; // Initial sum (without carry)
    // Propagate carry without creating additional Morton objects
    while (carry.any())
    {
        Morton<Dim, N> shiftedCarry = carry;
        shiftedCarry <<= 1;
        Morton<Dim, N> newCarry = result & shiftedCarry;
        result ^= shiftedCarry;
        carry = newCarry;
    }
    return result;
}

// Subtraction operators
template <size_t Dim, size_t N>
Morton<Dim,N> operator-(const Morton<Dim, N> &lhs, const Morton<Dim, N> &rhs)
{
    auto ones = static_cast<Morton<Dim, N>>((static_cast<size_t>(1) << Dim) - 1);
    auto complement = (~rhs + ones);
    return (lhs + complement);
}

/*
template <size_t Dim, size_t N>
Morton<Dim,N> operator<<(const Morton<Dim,N> &m, size_t shift)
{
    Morton<Dim,N> result(m);
    return (result.operator<<=(shift));
}
*/

template <size_t Dim, size_t N>
class MortonEncoder
{
public:
    MortonEncoder() { computeMagicBitsEnc(); }
    template <typename T>
    Morton<Dim, N> encode(T coord, size_t dir) const;
    template <typename T>
    Morton<Dim, N> encode(const std::array<T, Dim> &coord) const
    {
        return encodeImpl(coord, std::make_index_sequence<Dim>{});
    }
    BitArray<N> decode(Morton<Dim, N> m, size_t dir) const;
    std::array<BitArray<N>, Dim> decode(const Morton<Dim, N> &m) const
    {
        return decodeImpl(m, std::make_index_sequence<Dim>{});
    }

private:
    void computeMagicBitsEnc();
    template <typename T, size_t... Is>
    auto encodeImpl(const std::array<T, Dim> &coord, std::index_sequence<Is...>) const
    {
        Morton<Dim, N> m = (encode(coord[Is], Is) | ...);
        return m;
    }
    template <size_t... Is>
    auto decodeImpl(const Morton<Dim, N> &m, std::index_sequence<Is...>) const
    {
        return std::array<BitArray<N>, Dim>{{decode(m, Is)...}};
    }
    static constexpr size_t maxMove_ = (Dim - 1) * N;
    static constexpr size_t determineNumShift()
    {
        auto m = maxMove_;
        size_t i(0);
        for (; m != 0; ++i)
            m = (m >> 1);
        return i;
    }
    static constexpr size_t numShifts_ = determineNumShift();
    static constexpr size_t numMagicBits_ = numShifts_ + 1;
    static constexpr array<size_t, numMagicBits_> createShifts()
    {
        array<size_t, numMagicBits_> shfts;
        shfts[0] = (1 << numShifts_);
        for (auto i(1); i < numMagicBits_; ++i)
            shfts[i] = (shfts[i - 1] >> 1);
        return shfts;
    }
    static constexpr array<size_t, numMagicBits_> shifts_ = createShifts();
    array<BitArray<Dim*N>, numMagicBits_> magicBits_;
};

template <size_t Dim, size_t N>
void MortonEncoder<Dim, N>::computeMagicBitsEnc()
{
    for (auto i(0); i < N; ++i)
    {
        magicBits_[0].set(i);
    }
    array<BitArray<numShifts_>, N> shiftBits;
    for (auto i(0); i < N; ++i)
    {
        shiftBits[i] = BitArray<numShifts_>(i * (Dim - 1));
    }
    for (auto i(0); i < numShifts_; ++i)
    {
        for (auto j(0); j < N; ++j)
        {
            magicBits_[numShifts_ - i][j] = shiftBits[j][i];
        }
        magicBits_[numShifts_ - i] &= magicBits_[0];
    }
    auto maskLoc = magicBits_[0];
    unsigned int shft = (1 << (numShifts_ - 1));
    for (auto i(1); i < numMagicBits_; ++i)
    {
        for (auto j(i + 1); j < numMagicBits_; ++j)
            magicBits_[j] = (magicBits_[j] & ~magicBits_[i]) | ((magicBits_[j] & magicBits_[i]) << shft);
        maskLoc = (maskLoc & ~magicBits_[i]) | ((maskLoc & magicBits_[i]) << shft);
        magicBits_[i] = ((~magicBits_[i] & maskLoc) | (magicBits_[i] << shft & maskLoc));
        shft = (shft >> 1);
    }
};

template <size_t Dim, size_t N>
template <typename T>
Morton<Dim, N> MortonEncoder<Dim, N>::encode(T coord, size_t dir) const
{
    Morton<Dim,N> m(coord);
    BitArray<Dim*N>& b = reinterpret_cast<BitArray<Dim*N>&>(m);
    b &= magicBits_[0];
    for (auto i(1); i < numMagicBits_; ++i)
    {
        b = (b | (b << shifts_[i])) & magicBits_[i];
    }
    if (dir != 0)
    {
        b <<= (dir);
    }
    return m;
};

template <size_t Dim, size_t N>
BitArray<N> MortonEncoder<Dim, N>::decode(Morton<Dim, N> m, size_t dir) const
{
    BitArray<Dim*N>& b = reinterpret_cast<BitArray<Dim*N>&>(m);
    if (dir != 0)
    {
        b >>= dir;
    }
    b &= magicBits_[numMagicBits_ - 1];
    for (auto i(numMagicBits_ - 2); i != -1; --i)
    {
        b = ((b>> shifts_[i + 1]) | b) & (magicBits_[i]);
    }
    return m;
}

#endif // MORTON_HPP