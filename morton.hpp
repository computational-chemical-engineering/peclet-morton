#ifndef MORTON_HPP
#define MORTON_HPP 

#include <array>
#include <limits>
#include <cstring>
#include <utility>
#include "bitarray.hpp"

//using BitArray;
using std::array;

template<size_t Dim, size_t N>
class Morton : public BitArrayBase<Morton<Dim,N>, Dim*N> {
public:
    // Bring constructors from BitArrayBase into BitArray
    using BitArrayBase<Morton<Dim,N>, Dim*N>::BitArrayBase;

    // Constructor to convert from BitArrayBase to Morton
    Morton(const BitArrayBase<Morton<Dim, N>, Dim*N>& other)
        : BitArrayBase<Morton<Dim, N>, Dim*N>(other) {}

    Morton<Dim, N>& operator<<=(size_t pos)
    {
        BitArrayBase<Morton<Dim,N>, Dim*N>::operator<<= (Dim*pos);
        return *this;
    }

    template<typename T>
    Morton<Dim, N>& operator=(T val) {
        BitArrayBase<Morton<Dim,N>, Dim*N>::operator= (val);
        return *this;
    }

    Morton<Dim, N>& operator>>=(size_t pos)
    {
        BitArrayBase<Morton<Dim,N>, Dim*N>::operator>>= (Dim*pos);
        return *this;
    }

    // Addition operators
    Morton operator+(const Morton &other) const
    {
        Morton<Dim, N> carry = (*this) & other;  // Initial carry
        Morton<Dim, N> result = (*this) ^ other; // Initial sum (without carry)
        // Propagate carry
        while (carry.any())
        {
            Morton<Dim, N> shiftedCarry = (carry << 1);
            carry = result & shiftedCarry;
            result ^= shiftedCarry;
        }
        return result;
    }

    Morton& operator+=(const Morton& other) {
        Morton<Dim, N> carry = (*this) & other;  // Initial carry
        (*this) ^= other; // Initial sum (without carry)
        // Propagate carry
        while (carry.any())
        {
            Morton<Dim, N> shiftedCarry = (carry << 1);
            carry = (*this) & shiftedCarry;
            (*this) ^= shiftedCarry;
        }
        return *this;
    }

    // Subtraction operators
    Morton operator-(const Morton& other) const {
        auto ones = static_cast<Morton<Dim,N> >((static_cast<size_t>(1) << Dim)-1);
        auto complement = (~other + ones);
        return operator+(complement);
    }

    Morton& operator-=(const Morton& other) {
        Morton<Dim,N> ones = static_cast<Morton<Dim,N> >((static_cast<size_t>(1) << Dim)-1);
        Morton<Dim,N> complement = (~other + ones);
        return operator+=(complement);
    }

    bool test( size_t dir, size_t pos) const
    {
        return BitArrayBase<Morton<Dim,N>, Dim*N>::test(Dim*pos + dir);
    }

    Morton& set( size_t dir, size_t pos, bool value = true)
    {
        BitArrayBase<Morton<Dim,N>, Dim*N>::set(Dim*pos + dir, value);
        return *this;
    }   

    Morton& reset( size_t dir, size_t pos)
    {
        BitArrayBase<Morton<Dim,N>, Dim*N>::reset(Dim*pos + dir);
        return *this;
    }   

    static constexpr size_t szInBytes = (Dim*N+7)/8;
    static constexpr size_t szBitsetInBytes = (N+7)/8;

};

template <size_t Dim, size_t N>
Morton<Dim, N> operator<<(const Morton<Dim, N>& m, size_t shift) {
    Morton<Dim, N> result = m; // Create a copy of the input Morton object
    result <<= shift;          // Apply the left shift
    return result;             // Return the modified copy
}

template <size_t Dim, size_t N>
Morton<Dim, N> operator>>(const Morton<Dim, N>& m, size_t shift) {
    Morton<Dim, N> result = m; // Create a copy of the input Morton object
    result >>= shift;          // Apply the left shift
    return result;             // Return the modified copy
}


/*
template <size_t Dim, size_t N>
Morton<Dim,N> operator<<(const Morton<Dim,N> &m, size_t shift)
{
    Morton<Dim,N> result(m);
    return (result.operator<<=(shift));
}
*/

template<size_t Dim, size_t N>
class MortonEncoder {
public:
    MortonEncoder() {computeMagicBitsEnc();}
    template<size_t dir, typename T>
    Morton<Dim, N> encode(T coord) const;
    template<typename T>
    Morton<Dim, N> encode(const std::array<T, Dim>& coord) const {
        return encodeImpl(coord, std::make_index_sequence<Dim>{});
    }
    template<size_t dir> 
    BitArray<N> decode(Morton<Dim, N> m) const;
    std::array<BitArray<N>, Dim> decode(const Morton<Dim, N>& m) {
        return decodeImpl(m, std::make_index_sequence<Dim>{});
    }
private:
    void computeMagicBitsEnc();
    template <typename T, size_t... Is>
    auto encodeImpl(const std::array<T, Dim>& coord, std::index_sequence<Is...>) const {
        Morton<Dim, N> m = (encode<Is>(coord[Is]) | ...);
        return m;
    }
    template <size_t... Is>
    auto decodeImpl(const Morton<Dim, N>& m, std::index_sequence<Is...>) const {
        return std::array<BitArray<N>, Dim>{{ decode<Is>(m)... }};
    }
    static constexpr size_t maxMove_ = (Dim-1)*N;
    static constexpr size_t determineNumShift()
    {
            auto m = maxMove_;
            size_t i(0);
            for(; m !=0; ++i)
                m = (m >>1);
            return i;
    }
    static constexpr size_t numShifts_ = determineNumShift();
    static constexpr size_t numMagicBits_ = numShifts_ + 1;
    static constexpr array<size_t, numMagicBits_> createShifts() 
    {
        array<size_t, numMagicBits_> shfts;
        shfts[0] = (1<<numShifts_);
        for(auto i(1);i < numMagicBits_; ++i)
            shfts[i] = (shfts[i-1] >> 1);
        return shfts;
    }
    static constexpr array<size_t, numMagicBits_> shifts_ = createShifts();
    array<Morton<Dim,N>, numMagicBits_> magicBits_;
};
 
template<size_t Dim, size_t N>
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

template<size_t Dim, size_t N>
template<size_t dir, typename T> 
Morton<Dim, N> MortonEncoder<Dim, N>::encode(T coord) const
{
    static_assert(dir<Dim, "direction should be in between 0 and Dim.");
    Morton<Dim, N> m(coord);
    m &= magicBits_[0];
    for(auto i(1); i < numMagicBits_; ++i){
        Morton<Dim, N> m2(m);
        m2.BitArrayBase<Morton<Dim,N>, Dim*N>::operator<<=(shifts_[i]);
        m.BitArrayBase<Morton<Dim,N>, Dim*N>::operator|=(m2);
        m.BitArrayBase<Morton<Dim,N>, Dim*N>::operator&=(magicBits_[i]);
    }
    if constexpr(dir != 0) {
        m.BitArrayBase<Morton<Dim,N>, Dim*N>::operator<<=(dir);
    }
    return m;
         //     bs = (bs | (bs << shifts_[i])) & magicBits_[i];
};

template<size_t Dim, size_t N>
template<size_t dir>
BitArray<N> MortonEncoder<Dim, N>::decode(Morton<Dim, N> m) const
{
    static_assert(dir<Dim, "direction should be in between 0 and Dim.");
    if constexpr(dir != 0) {
        m.BitArrayBase<Morton<Dim,N>, Dim*N>::operator>>=(dir);
    }
    m &= magicBits_[numMagicBits_-1];
    for(auto i(numMagicBits_-2); i != -1; --i){
        Morton<Dim, N> m2(m);
        m2.BitArrayBase<Morton<Dim,N>, Dim*N>::operator>>(shifts_[i+1]);
        m.BitArrayBase<Morton<Dim,N>, Dim*N>::operator|=(m2);
        m.BitArrayBase<Morton<Dim,N>, Dim*N>::operator&=(magicBits_[i]);
    }
    BitArray<N> b(m);
    //std::memcpy(&b, &m, Morton<Dim, N>::szBitsetInBytes);
    return b;
}

#endif //MORTON_HPP