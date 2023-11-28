#ifndef MORTON_HPP
#define MORTON_HPP 

#include <bitset>
#include <array>
#include <limits>
#include <cstring>

using std::bitset;
using std::array;

template<size_t Dim, size_t N>
class Morton: public bitset<Dim*N> {
public:
    template<typename T>
    constexpr Morton(T val) : std::bitset<Dim*N>(val) {}

    Morton operator<<(size_t pos) const
    {
        return Morton(this->std::bitset<Dim*N>::operator<<(Dim * pos));
    }

    Morton<Dim, N>& operator<<=(size_t pos)
    {
        this->bitset<Dim*N>::operator<<= (Dim*pos);
        return *this;
    }

    template<typename T>
    Morton<Dim, N>& operator=(T val) {
        this->bitset<Dim*N>::operator= (val);
        return *this;
    }

    Morton operator>>(size_t pos) const
    {
        return Morton(this->std::bitset<Dim*N>::operator>>(Dim * pos));
    }

    Morton<Dim, N>& operator>>=(size_t pos)
    {
        this->bitset<Dim*N>::operator>>= (Dim*pos);
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
            Morton<Dim, N> shiftedCarry = carry << 1;
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
            Morton<Dim, N> shiftedCarry = carry << 1;
            carry = (*this) & shiftedCarry;
            (*this) ^= shiftedCarry;
        }
        return *this;
    }

    // Subtraction operators
    Morton operator-(const Morton& other) const {
        Morton<Dim, N> complement = ~other;
        complement += ones;
        return (*this) + complement;
    }

    Morton& operator-=(const Morton& other) {
        Morton<Dim, N> complement = ~other;
        complement += ones;
        (*this) += complement;
        return *this;
    }

    bool operator<(const Morton& other) const;

    bool operator>(const Morton& other) const{
        return other.operator>(*this);
    }

    bool operator>=(const Morton& other) const{
        return !(this->operator<(*this));
    }

    bool operator<=(const Morton& other) const{
        return !(this->operator>(*this));
    }

    using bitset<Dim*N>::test;

    bool test( size_t dir, size_t pos) const
    {
        return this->test(Dim*pos + dir);
    }

    Morton& set( size_t dir, size_t pos)
    {
        this->bitset<Dim*N>::set(Dim*pos + dir, true);
        return *this;
    }   

    Morton& reset( size_t dir, size_t pos)
    {
        this->bitset<Dim*N>::reset(Dim*pos + dir);
        return *this;
    }   

    Morton& flip( size_t dir, size_t pos)
    {
        this->bitset<Dim*N>::flip(Dim*pos + dir);
        return *this;
    }

    static constexpr size_t szInBytes = (Dim*N+7)/8;
    static constexpr size_t szBitsetInBytes = (N+7)/8;

    private:
    static constexpr Morton<Dim, N> createOnes() {
        size_t i(1);
        i <<= Dim;
        i-=1;
        return static_cast<Morton<Dim,N> >(i);
    }
    public:
    static constexpr Morton<Dim,N> ones = createOnes();
};
 
template<size_t Dim, size_t N>
class MortonEncoder {
public:
    MortonEncoder() {computeMagicBitsEnc();}
    template<size_t dir, typename T> 
    Morton<Dim, N> encode(T coord) const;
    template<typename T>
    Morton<Dim, N> encode(array<T,Dim> coord) const{
        Morton<Dim, N> m(encode<0>(coord[0]));
        for(auto i(1); i<Dim;++i)
            m |= encode<i>(coord[0]);
        return m;
    }
    template<size_t dir> 
    bitset<N> decode(Morton<Dim, N> m) const;
    array<bitset<N>,Dim> decode(Morton<Dim, N> m) const
    {
        array<bitset<N>,Dim> coord;
        for(auto i(0); i<Dim;++i)
            coord[i] = decode<i>(m);
        return coord;
    }
private:
    void computeMagicBitsEnc();
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
    array<bitset<Dim*N>, numMagicBits_> magicBits_;
};
 
template<size_t Dim, size_t N>
void MortonEncoder<Dim, N>::computeMagicBitsEnc()
{
    for (auto i(0); i < N; ++i)
    {
        magicBits_[0].set(i);
    }
    array<bitset<numShifts_>, N> shiftBits;
    for (auto i(0); i < N; ++i)
    {
        shiftBits[i] = bitset<numShifts_>(i * (Dim - 1));
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
        Morton<Dim, N> m2(m.bitset<Dim*N>::operator<<(shifts_[i]));
        m.bitset<Dim*N>::operator|=(m2);
        m.bitset<Dim*N>::operator&=(magicBits_[i]);
    }
    if constexpr(dir != 0) {
        m.bitset<Dim*N>::operator<<=(dir);
    }
    return m;
         //     bs = (bs | (bs << shifts_[i])) & magicBits_[i];
};

template<size_t Dim, size_t N>
template<size_t dir>
bitset<N> MortonEncoder<Dim, N>::decode(Morton<Dim, N> m) const
{
    static_assert(dir<Dim, "direction should be in between 0 and Dim.");
    if constexpr(dir != 0) {
        m.bitset<Dim*N>::operator>>=(dir);
    }
    m &= magicBits_[numMagicBits_-1];
    for(auto i(numMagicBits_-2); i != -1; --i){
        Morton<Dim, N> m2(m.bitset<Dim*N>::operator>>(shifts_[i+1]));
        m.bitset<Dim*N>::operator|=(m2);
        m.bitset<Dim*N>::operator&=(magicBits_[i]);
    }
    bitset<N> b;
    std::memcpy(&b, &m, Morton<Dim, N>::szBitsetInBytes);
    return b;
}

template<size_t Dim, size_t N>
bool Morton<Dim, N>::operator<(const Morton& other) const
{
    constexpr size_t numBits = (Dim*N);
    if constexpr (numBits <= sizeof(unsigned long) * 8)
    {
        // If the bitset fits into an unsigned long, use to_ulong for comparison
        return this->to_ulong() < other.to_ulong();
    }
    else if constexpr (numBits <= sizeof(unsigned long long) * 8)
    {
        // If the bitset fits into an unsigned long long, use to_ullong for comparison
        return this->to_ullong() < other.to_ullong();
    }
    else
    {
        constexpr unsigned long max_ulong = std::numeric_limits<unsigned long>::max();
        constexpr std::bitset<numBits> mask_ulong(max_ulong);
        constexpr size_t num_ulong_bits = sizeof(unsigned long) * 8;
        constexpr size_t iterations = (numBits + num_ulong_bits - 1) / num_ulong_bits; // Calculate number of segments
        unsigned long segments_a[iterations], segments_b[iterations];
        std::bitset<numBits> temp_a = *this, temp_b = other;
        for (size_t i = 0; i < iterations; ++i)
        {
            segments_a[i] = (temp_a & mask_ulong).to_ulong();
            temp_a >>= num_ulong_bits;
            segments_b[i] = (temp_b & mask_ulong).to_ulong();
            temp_b >>= num_ulong_bits;
        }
        for (size_t i = iterations - 1; i != (-1); --i)
            if (segments_a[i] != segments_b[i])
            {
                return segments_a[i] < segments_b[i];
            }
        return false; // Equal if all segments are the same
    }
};

// Custom hash functor for Morton
namespace std {
    template <size_t Dim, size_t N>
    struct hash<Morton<Dim, N>> {
        size_t operator()(const Morton<Dim, N>& m) const {
            return std::hash<std::bitset<Dim * N>>{}(reinterpret_cast<const std::bitset<Dim * N>&>(m));
        }
    };
}

#endif //MORTON_HPP