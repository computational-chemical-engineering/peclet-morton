#ifndef BITARRAY_HPP
#define BITARRAY_HPP

#include <cstring> // For memcpy
#include <cstdint>
#include <array>
#include <type_traits>
#include <algorithm>
#include <stdexcept>

template <typename Derived, size_t N>
class BitReference;

template <typename Derived, size_t N>
class BitArrayBase
{
public:
    constexpr BitArrayBase()
    {
//        std::memset(bits_, 0, sizeof(bits_));
        std::fill(std::begin(bits_), std::end(bits_), static_cast<IntType>(0));
    }
    // Copy constructor
    constexpr BitArrayBase(const Derived &other)
    {
//        std::memcpy(bits_, other.bits_, sizeof(bits_));
        std::copy(std::begin(other.bits_), std::end(other.bits_), std::begin(bits_));
    }

    // Template constructor to initialize from BitArrayBase<M>
    template <typename OtherDerived, size_t M>
    BitArrayBase(const BitArrayBase<OtherDerived, M>& other)
    {
        std::memset(bits_, 0, sizeof(bits_));
        constexpr size_t minSizeBytes = (std::min(N, M) + 7) / 8;
        std::memcpy(bits_, other.bits_, minSizeBytes);
    }

    template <typename T>
    BitArrayBase(T value)
    {
        std::memset(bits_, 0, sizeof(bits_));
        constexpr size_t minSizeBytes = std::min(sizeof(bits_), sizeof(value));
        std::memcpy(bits_, &value, minSizeBytes);
    }

    Derived &operator=(const Derived &other)
    {
        if (this != &other)
        { // Check for self-assignment
            std::copy(std::begin(other.bits_), std::end(other.bits_), std::begin(bits_));
        }
        return *this;
    }

    bool test(size_t index) const
    {
        if (index >= N)
        {
            throw std::out_of_range("BitArrayBase index out of range");
        }
        if constexpr (numInt == 1)
        {
            return (bits_[0] & (static_cast<IntType>(1) << index)) != 0;
        }
        else
        {
            size_t arrayIndex = index / szInt;  // Calculate which integer in the array
            size_t bitPosition = index % szInt; // Calculate the position within that integer
            return (bits_[arrayIndex] & (static_cast<IntType>(1) << bitPosition)) != 0;
        }
    }

    BitReference<Derived, N> operator[](size_t index) {
        return BitReference<Derived, N>(*this, index);
    }

    bool operator[](size_t index) const
    {
        return test(index);
    }

    void set(size_t index, bool value=true)
    {
        if (index >= N)
        {
            throw std::out_of_range("BitArrayBase index out of range");
        }
        if constexpr (numInt == 1)
        {
            bits_[0] |= static_cast<IntType>(value) << index;
        }
        else
        {
            size_t arrayIndex = index / szInt;  // Calculate which integer in the array
            size_t bitPosition = index % szInt; // Calculate the position within that integer
            bits_[arrayIndex] |= static_cast<IntType>(value) << bitPosition;
        }
    }

    void reset(size_t index)
    {
        if (index >= N)
        {
            throw std::out_of_range("BitArrayBase index out of range");
        }
        if constexpr (numInt == 1)
        {
            bits_[0] &= ~(static_cast<IntType>(1) << index);
        }
        else
        {
            size_t arrayIndex = index / szInt;  // Calculate which integer in the array
            size_t bitPosition = index % szInt; // Calculate the position within that integer
            bits_[arrayIndex] &= ~(static_cast<IntType>(1)) << bitPosition;
        }
    }

/*
    // Friend declaration for operator<<
    template <typename Derived, size_t N>
    friend std::ostream &operator<<(std::ostream &os, const BitArrayBase<Derived,N> &bitArray)
    {
        for (size_t i = N; i > 0; --i)
        {
            os << bitArray.test(i - 1); // Output each bit, test() method is assumed to be implemented
        }
        return os;
    }
    // Right shift operator
    Derived &operator>>=(size_t shift)
    {
        if (shift >= N)
        {
            std::fill(std::begin(bits_), std::end(bits_), 0); // Clear all bits_ if shift >= N
            return *this;
        }
        if constexpr (numInt == 1)
        {
            bits_[0] >>= shift;
        }
        else
        {
            size_t wholeShifts = shift / szInt;
            size_t bitShifts = shift % szInt;

            // Shift whole integers
            if (wholeShifts != 0)
            {
                for (size_t i = 0; i < (numInt - wholeShifts); ++i)
                {
                    bits_[i] = bits_[i + wholeShifts];
                }
                // Zero out the vacated integers
                for (size_t i = (numInt - wholeShifts); i < numInt; ++i)
                {
                    bits_[i] = 0;
                }
            }
            for (size_t i = 0; i < numInt - 1; ++i)
            {
                bits_[i] = (bits_[i] >> bitShifts) | (bits_[i + 1] << (szInt - bitShifts));
            }
            bits_[numInt - 1] >>= bitShifts; // Shift the last integer
        }
        return *this;
    }
*/
    Derived &operator<<=(size_t shift)
    {
        if (shift >= N)
        {
            std::fill(std::begin(bits_), std::end(bits_), 0); // Clear all bits_ if shift >= N
            return *this;
        }
        if constexpr (numInt == 1)
        {
            bits_[0] <<= shift;
            // Mask out bits_ that are outside the range of N
            if constexpr (N < szInt)
            {
                bits_[0] &= ((static_cast<IntType>(1) << N) - 1);
            }
        }
        else
        {
            size_t wholeShifts = shift / szInt;
            size_t bitShifts = shift % szInt;

            // Shift whole integers
            if (wholeShifts != 0)
            {
                for (size_t i = numInt - 1; i >= wholeShifts; --i)
                {
                    bits_[i] = bits_[i - wholeShifts];
                }
                // Zero out the vacated integers
                for (size_t i = 0; i < wholeShifts; ++i)
                {
                    bits_[i] = 0;
                }
            }

            // Shift individual bits_
            if (bitShifts > 0)
            {
                for (size_t i = numInt - 1; i > 0; --i)
                {
                    bits_[i] = (bits_[i] << bitShifts) | (bits_[i - 1] >> (szInt - bitShifts));
                }
                bits_[0] <<= bitShifts;
            }

            // Mask out any excess bits_ in the last storage element
            if constexpr (N % szInt != 0)
            {
                bits_[numInt - 1] &= ((static_cast<IntType>(1) << (N % szInt)) - 1);
            }
        }
        return *this;
    }

    Derived &operator|=(const Derived &other)
    {
        if constexpr (numInt == 1)
        {
            bits_[0] |= other.bits_[0];
        }
        else
        {
            for (size_t i = 0; i < numInt; ++i)
            {
                bits_[i] |= other.bits_[i];
            }
        }
        return *this;
    }

    // Bitwise AND assignment operator
    Derived &operator&=(const Derived &other)
    {
        if constexpr (numInt == 1)
        {
            bits_[0] &= other.bits_[0];
        }
        else
        {
            for (size_t i = 0; i < numInt; ++i)
            {
                bits_[i] &= other.bits_[i];
            }
        }
        return *this;
    }

    // Bitwise XOR assignment operator
    Derived &operator^=(const Derived &other)
    {
        if constexpr (numInt == 1)
        {
            bits_[0] ^= other.bits_[0];
        }
        else
        {
            for (size_t i = 0; i < numInt; ++i)
            {
                bits_[i] ^= other.bits_[i];
            }
        }
        return *this;
    }

    // Unary NOT operator
    Derived operator~() const
    {
        Derived result;

        if constexpr (numInt == 1)
        {
            result.bits_[0] = ~bits_[0];
        }
        else
        {
            for (size_t i = 0; i < numInt; ++i)
            {
                result.bits_[i] = ~bits_[i];
            }
        }
        // If N is not a multiple of the size of IntType in bits_, mask off the excess bits_
        constexpr size_t excessBits = N % (sizeof(IntType) * 8);
        if (excessBits != 0)
        {
            result.bits_[numInt - 1] &= (static_cast<IntType>(1) << excessBits) - 1;
        }
        return result;
    }

    bool operator<(const Derived& other) const
    {
        if constexpr (numInt == 1)
        {
            return (bits_[0] < other.bits_[0]);
        } else {
            for (size_t i = numInt; i >0 ; --i)
            {
                if (bits_[i-1] < other.bits_[i-1])
                    return true;
            }
        }
        return false;
    }
    bool operator>=(const Derived& other) const
    {
        return ~(operator<(other));
    }
    bool operator>(const Derived& other) const
    {
        return (other.operator<(*this));
    }
    bool operator<=(const Derived& other) const
    {
        return ~(other.operator>(*this));
    }

    bool operator==(const Derived& other) const
    {
        if constexpr (numInt == 1)
        {
            return (bits_[0] == other.bits_[0]);
        } else {
            for (size_t i = 0; i < numInt ; ++i)
            {
                if (bits_[i] != other.bits_[i])
                    return false;
            }
        }
        return true;
    }
    bool operator!=(const Derived& other) const
    {
        return ~(operator==(other));
    }
private:
    static constexpr size_t szInt =
        N <= 8 ? 8 : N <= 16 ? 16
                 : N <= 32   ? 32
                 : N <= 64   ? 64
                 : N <= 96   ? 32
                             : 64;

    static constexpr size_t numInt = (N + szInt - 1) / szInt;

    using IntType =
        std::conditional_t<szInt == 8, uint8_t,
                           std::conditional_t<szInt == 16, uint16_t,
                                              std::conditional_t<szInt == 32, uint32_t,
                                                                 uint64_t>>>; // Default to uint64_t

    IntType bits_[numInt];
};

template <size_t N>
class BitArray : public BitArrayBase<BitArray<N>, N> {
public:
    // Bring constructors from BitArrayBase into BitArray
    using BitArrayBase<BitArray<N>, N>::BitArrayBase;
};

template <typename Derived, size_t N>
class BitReference {
    BitArrayBase<Derived, N>& bitArray;
    size_t index;

public:
    BitReference(BitArrayBase<Derived, N>& ba, size_t idx) : bitArray(ba), index(idx) {}

    // Assignment operator to set the bit
    BitReference& operator=(bool value) {
        // Assuming Derived implements a set method
        static_cast<Derived&>(bitArray).set(index, value);
        return *this;
    }

    // Conversion to bool (read the bit)
    operator bool() const {
        // Assuming Derived implements a test method
        return static_cast<const Derived&>(bitArray).test(index);
    }
};


template <typename Derived, size_t N>
BitArrayBase<Derived, N> operator>>(const BitArrayBase<Derived, N> &bitArray, size_t shift)
{
    BitArrayBase<Derived, N> result = bitArray; // Initialize result with the value of bitArray
    result >>= shift;              // Apply right shift using operator>>=
    return result;                 // Return the result
}

template <typename Derived, size_t N>
BitArrayBase<Derived, N> operator<<(const BitArrayBase<Derived, N> &bitArray, size_t shift)
{
    BitArrayBase<Derived, N> result = bitArray; // Initialize result with the value of bitArray
    result <<= shift;              // Apply right shift using operator>>=
    return result;                 // Return the result
}

template <typename Derived, size_t N>
Derived operator|(const BitArrayBase<Derived, N> &bitArray, const BitArrayBase<Derived, N> &bitArray2)
{
    auto result = bitArray;
    result |= bitArray2;
    return result;
}
template <typename Derived, size_t N>
Derived operator&(const BitArrayBase<Derived, N> &bitArray, const BitArrayBase<Derived, N> &bitArray2)
{
    auto result = bitArray;
    result &= bitArray2;
    return result;
}
template <typename Derived, size_t N>
Derived operator^(const BitArrayBase<Derived, N> &bitArray, const BitArrayBase<Derived, N> &bitArray2)
{
    auto result = bitArray;
    result ^= bitArray2;
    return result;
}
#endif // BITARRAY_HPP