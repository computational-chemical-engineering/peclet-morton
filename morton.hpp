#ifndef MORTON_HPP
#define MORTON_HPP 

#include <bitset>
#include <limits>

using std::bitset;

template<size_t Dim, size_t N>
class Morton: public bitset<Dim*N> {
public:
    using bitset<Dim*N>::bitset;

    Morton(const bitset<Dim*N>& bs) : bitset<Dim*N>(bs) {}

    Morton operator<<(size_t pos) const
    {
        return Morton(this->std::bitset<Dim*N>::operator<<(Dim * pos));
    }

    Morton<Dim, N>& operator<<=(size_t pos)
    {
        this->bitset<Dim*N>::operator<<= (Dim*pos);
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

    bool operator<(const Morton& other);

    bool operator>(const Morton& other){
        return other.operator>(*this);
    }

    bool operator>=(const Morton& other){
        return !(this->operator<(*this));
    }

    bool operator<=(const Morton& other){
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
bool Morton<Dim, N>::operator<(const Morton& other)
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
            return std::hash<std::bitset<Dim * N>>{}(static_cast<const std::bitset<Dim * N>&>(m));
        }
    };
}

#endif //MORTON_HPP 

