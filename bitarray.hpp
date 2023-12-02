#ifndef BITARRAY_HPP
#define BITARRAY_HPP

#include <cstring>
#include <cstdint>
#include <array>
#include <type_traits>
#include <algorithm>
#include <stdexcept>
#include <bit>

// Forward declaration for BitReference
template <typename Derived, size_t N>
class BitReference;

/**
 * @brief A base class template for creating bit array types.
 * 
 * This class template provides a flexible base for implementing bit arrays of a fixed size.
 * It supports various bit manipulation operations, bitwise operators, comparison operators,
 * and also includes support for streaming to standard output.
 * It is very similar to std::bitset, but this was created because of the need of more control
 * in the derived class Morton.
 * 
 * @tparam Derived The derived type that inherits from BitArrayBase (CRTP Pattern).
 * @tparam N The size of the bit array in bits.
 */
template <typename Derived, size_t N>
class BitArrayBase
{
public:
    /**
     * @brief Default constructor.
     * Initializes all bits to 0.
     */
    constexpr BitArrayBase();

    /**
     * @brief Copy constructor.
     * 
     * @param other The BitArrayBase object to copy from.
     */
    constexpr BitArrayBase(const Derived &other);

    /**
     * @brief Template constructor to initialize from another BitArrayBase of potentially different size.
     * 
     * @tparam OtherDerived The derived type of the other BitArrayBase.
     * @tparam M The size of the other BitArrayBase.
     * @param other The other BitArrayBase object to initialize from.
     */
    template <typename OtherDerived, size_t M>
    BitArrayBase(const BitArrayBase<OtherDerived, M> &other);

    /**
     * @brief Constructor to initialize the bit array from an integer value.
     * 
     * @tparam T The integer type.
     * @param value The integer value to initialize the bit array with.
     */
    template <typename T>
    BitArrayBase(T value);

    /**
     * @brief Assignment operator.
     * 
     * @param other The BitArrayBase object to assign from.
     * @return A reference to the modified object.
     */
    Derived &operator=(const Derived &other);

 /**
 * @brief Tests a specific bit for its value.
 * 
 * @param index The index of the bit to test.
 * @return True if the bit is set (1), false otherwise.
 * @exception std::out_of_range If the index is out of range.
 */
bool test(size_t index) const;

/**
 * @brief Access operator for const objects.
 * 
 * @param index The index of the bit to access.
 * @return The value of the bit at the given index.
 */
bool operator[](size_t index) const;

/**
 * @brief Return the bitarray as an integer value
 *
 * @tparam T The integer type.
 * @return value The integer value to initialize the bit array with.
 */
template <typename T>
T getIntValue() const;

/**
 * @brief Sets a specific bit to 1.
 * 
 * @param index The index of the bit to set.
 * @exception std::out_of_range If the index is out of range.
 */
void set(size_t index);

/**
 * @brief Sets a specific bit to a given value.
 * 
 * @param index The index of the bit to set.
 * @param value The value to set the bit to (true for 1, false for 0).
 * @exception std::out_of_range If the index is out of range.
 */
void set(size_t index, bool value);

/**
 * @brief Resets a specific bit to 0.
 * 
 * @param index The index of the bit to reset.
 * @exception std::out_of_range If the index is out of range.
 */
void reset(size_t index);

/**
 * @brief Right shift assignment operator.
 * 
 * @param shift The number of positions to shift the bits to the right.
 * @return A reference to the modified object.
 */
Derived &operator>>=(size_t shift);

/**
 * @brief Left shift assignment operator.
 * 
 * @param shift The number of positions to shift the bits to the left.
 * @return A reference to the modified object.
 */
Derived &operator<<=(size_t shift);

/**
 * @brief Bitwise OR assignment operator.
 * 
 * @param other The BitArrayBase object to perform bitwise OR with.
 * @return A reference to the modified object.
 */
Derived &operator|=(const Derived &other);

/**
 * @brief Bitwise AND assignment operator.
 * 
 * @param other The BitArrayBase object to perform bitwise AND with.
 * @return A reference to the modified object.
 */
Derived &operator&=(const Derived &other);

/**
 * @brief Bitwise XOR assignment operator.
 * 
 * @param other The BitArrayBase object to perform bitwise XOR with.
 * @return A reference to the modified object.
 */
Derived &operator^=(const Derived &other);

/**
 * @brief Unary NOT operator.
 * 
 * @return A new object with all bits inverted.
 */
Derived operator~() const;

/**
 * @brief Less than comparison operator.
 * 
 * @param other The BitArrayBase object to compare with.
 * @return True if this object is less than the other, false otherwise.
 */
bool operator<(const Derived &other) const;

/**
 * @brief Greater than or equal comparison operator.
 * 
 * @param other The BitArrayBase object to compare with.
 * @return True if this object is greater than or equal to the other, false otherwise.
 */
bool operator>=(const Derived &other) const;

/**
 * @brief Greater than comparison operator.
 * 
 * @param other The BitArrayBase object to compare with.
 * @return True if this object is greater than the other, false otherwise.
 */
bool operator>(const Derived &other) const;

/**
 * @brief Less than or equal comparison operator.
 * 
 * @param other The BitArrayBase object to compare with.
 * @return True if this object is less than or equal to the other, false otherwise.
 */
bool operator<=(const Derived &other) const;

/**
 * @brief Equality comparison operator.
 * 
 * @param other The BitArrayBase object to compare with.
 * @return True if both objects are equal, false otherwise.
 */
bool operator==(const Derived &other) const;

/**
 * @brief Inequality comparison operator.
 * 
 * @param other The BitArrayBase object to compare with.
 * @return True if the objects are not equal, false otherwise.
 */
bool operator!=(const Derived &other) const;

/**
 * @brief Checks if all bits are set to 1.
 * 
 * @return True if all bits are set, false otherwise.
 */
bool all();

/**
 * @brief Checks if any bit is set to 1.
 * 
 * @return True if any bit is set, false otherwise.
 */
bool any();

/**
 * @brief Checks if no bits are set to 1.
 * 
 * @return True if no bits are set, false otherwise.
 */
bool none();

/**
 * @brief Returns the number of consecutive 0 bits starting from the least significant bit ("right").
 * 
 * @return Number of consecutive 0 bits 
 */
size_t countr_zero() const;

/**
 * @brief Provides a reference to a specific bit.
 * 
 * @param index The index of the bit to reference.
 * @return A BitReference object referring to the specified bit.
 */
BitReference<Derived, N> operator[](size_t index);

/**
 * @brief Stream output operator for BitArrayBase.
 * 
 * @tparam D The derived type.
 * @tparam M The size of the bit array.
 * @param os The output stream to write to.
 * @param bitArray The BitArrayBase object to output.
 * @return A reference to the output stream.
 */
template <typename D, size_t M>
friend std::ostream &operator<<(std::ostream &os, const BitArrayBase<D, M> &bitArray);

private:
    /**
     * @brief Size of the integer type used to store the bits, determined based on the size of the bit array.
     * 
     * It chooses the smallest integer size that can hold N bits, with a fallback to 64-bit integers.
     */
    static constexpr size_t szInt =
            N <= 8 ? 8 : N <= 16 ? 16
                 : N <= 32   ? 32
                 : N <= 64   ? 64
                 : N <= 96   ? 32
                             : 64;

    /**
     * @brief The number of integers of szInt size required to store N bits.
     * 
     * This value determines the size of the array used to store the bits.
     */
    static constexpr size_t numInt = (N + szInt - 1) / szInt;

    /**
     * @brief Alias for the underlying integer type used to store the bits.
     * 
     * The type is selected based on szInt to efficiently store the bits while minimizing memory usage.
     */
    using IntType =
        std::conditional_t<szInt == 8, uint8_t,
        std::conditional_t<szInt == 16, uint16_t,
        std::conditional_t<szInt == 32, uint32_t,
        uint64_t>>>; // Default to uint64_t if no smaller type can be used

    /**
     * @brief Array to store the bits of the bit array.
     * 
     * The array is of size numInt and each element is of type IntType.
     */
    IntType bits_[numInt];
};

/**
 * @brief Specialization of BitArrayBase for a specific size.
 * 
 * @tparam N The size of the bit array.
 */
template <size_t N>
class BitArray : public BitArrayBase<BitArray<N>, N>
{
public:
    /**
     * @brief Inherit constructors from BitArrayBase.
     */
    using BitArrayBase<BitArray<N>, N>::BitArrayBase;

    /**
     * @brief Construct a BitArray object from a BitArrayBase object.
     * 
     * @param other The BitArrayBase object to construct from.
     */
    BitArray(const BitArrayBase<BitArray<N>, N> &other);
};

/**
 * @brief A reference class to represent a single bit within a BitArray.
 * 
 * @tparam Derived The derived BitArray type.
 * @tparam N The size of the bit array.
 */
template <typename Derived, size_t N>
class BitReference
{
public:
    /**
     * @brief Construct a BitReference object.
     * 
     * @param ba Reference to the BitArrayBase object.
     * @param idx The index of the bit this BitReference refers to.
     */
    BitReference(BitArrayBase<Derived, N> &ba, size_t idx);

    /**
     * @brief Assignment operator to set the referenced bit.
     * 
     * @param value The boolean value to assign to the bit.
     * @return A reference to the BitReference object.
     */
    BitReference &operator=(bool value);

    /**
     * @brief Conversion operator to bool, to allow reading the value of the referenced bit.
     * 
     * @return The boolean value of the referenced bit.
     */
    operator bool() const;

private:
    /**
     * @brief Reference to the BitArrayBase object containing the bit.
     */
    BitArrayBase<Derived, N> &bitArray;

    /**
     * @brief The index of the bit within the BitArrayBase object.
     */
    size_t index;
};

/**
 * @brief Right shift operator for BitArray.
 * 
 * @tparam Derived The derived type of the BitArray.
 * @tparam N The size of the bit array.
 * @param bitArray The bit array to be shifted.
 * @param shift The number of positions to shift.
 * @return Derived The resulting bit array after the right shift operation.
 */
template <typename Derived, size_t N>
Derived operator>>(const BitArrayBase<Derived, N> &bitArray, size_t shift);

/**
 * @brief Left shift operator for BitArray.
 * 
 * @tparam Derived The derived type of the BitArray.
 * @tparam N The size of the bit array.
 * @param bitArray The bit array to be shifted.
 * @param shift The number of positions to shift.
 * @return Derived The resulting bit array after the left shift operation.
 */
template <typename Derived, size_t N>
Derived operator<<(const BitArrayBase<Derived, N> &bitArray, size_t shift);

/**
 * @brief Bitwise OR operator for BitArray.
 * 
 * @tparam Derived The derived type of the BitArray.
 * @tparam N The size of the bit array.
 * @param bitArray1 The first bit array operand.
 * @param bitArray2 The second bit array operand.
 * @return Derived The resulting bit array after the bitwise OR operation.
 */
template <typename Derived, size_t N>
Derived operator|(const BitArrayBase<Derived, N> &bitArray1, const BitArrayBase<Derived, N> &bitArray2);

/**
 * @brief Bitwise AND operator for BitArray.
 * 
 * @tparam Derived The derived type of the BitArray.
 * @tparam N The size of the bit array.
 * @param bitArray1 The first bit array operand.
 * @param bitArray2 The second bit array operand.
 * @return Derived The resulting bit array after the bitwise AND operation.
 */
template <typename Derived, size_t N>
Derived operator&(const BitArrayBase<Derived, N> &bitArray1, const BitArrayBase<Derived, N> &bitArray2);

/**
 * @brief Bitwise XOR operator for BitArray.
 * 
 * @tparam Derived The derived type of the BitArray.
 * @tparam N The size of the bit array.
 * @param bitArray1 The first bit array operand.
 * @param bitArray2 The second bit array operand.
 * @return Derived The resulting bit array after the bitwise XOR operation.
 */
template <typename Derived, size_t N>
Derived operator^(const BitArrayBase<Derived, N> &bitArray1, const BitArrayBase<Derived, N> &bitArray2);

// Implementations:

// BitArrayBase Default Constructor Implementation
template <typename Derived, size_t N>
constexpr BitArrayBase<Derived, N>::BitArrayBase():bits_{} {}

// BitArrayBase Copy Constructor Implementation
template <typename Derived, size_t N>
constexpr BitArrayBase<Derived, N>::BitArrayBase(const Derived &other)
{
    std::copy(std::begin(other.bits_), std::end(other.bits_), std::begin(bits_));
}

// BitArrayBase Template Constructor Implementation
template <typename Derived, size_t N>
template <typename OtherDerived, size_t M>
BitArrayBase<Derived, N>::BitArrayBase(const BitArrayBase<OtherDerived, M> &other):bits_{}
{
    constexpr size_t minSizeBytes = (std::min(N, M) + 7) / 8;
    std::memcpy(bits_, other.bits_, minSizeBytes);
}

// BitArrayBase Constructor from Value Implementation
template <typename Derived, size_t N>
template <typename T>
BitArrayBase<Derived, N>::BitArrayBase(T value):bits_{}
{
    constexpr size_t minSizeBytes = std::min(sizeof(bits_), sizeof(value));
    std::memcpy(bits_, &value, minSizeBytes);
}

// Assignment Operator Implementation
template <typename Derived, size_t N>
Derived &BitArrayBase<Derived, N>::operator=(const Derived &other)
{
    if (this != &other)
    { // Protect against self-assignment
        std::copy(std::begin(other.bits_), std::end(other.bits_), std::begin(bits_));
    }
    return *static_cast<Derived *>(this);
}

// Bit Test Method Implementation
template <typename Derived, size_t N>
bool BitArrayBase<Derived, N>::test(size_t index) const
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
        size_t arrayIndex = index / szInt;  // Determine which integer in the array
        size_t bitPosition = index % szInt; // Determine the bit position within that integer
        return (bits_[arrayIndex] & (static_cast<IntType>(1) << bitPosition)) != 0;
    }
}

// Operator[] Implementation
template <typename Derived, size_t N>
bool BitArrayBase<Derived, N>::operator[](size_t index) const
{
    return test(index);
}

template <typename Derived, size_t N>
template <typename T>
T BitArrayBase<Derived, N>::getIntValue() const
{
    T value;
    constexpr size_t minSizeBytes = std::min(sizeof(bits_), sizeof(value));
    std::memcpy(&value, bits_, minSizeBytes);
    return value;
}

// Set Bit Method Implementation
template <typename Derived, size_t N>
void BitArrayBase<Derived, N>::set(size_t index)
{
    if (index >= N)
    {
        throw std::out_of_range("BitArrayBase index out of range");
    }
    if constexpr (numInt == 1)
    {
        bits_[0] |= static_cast<IntType>(1) << index;
    }
    else
    {
        size_t arrayIndex = index / szInt;
        size_t bitPosition = index % szInt;
        bits_[arrayIndex] |= static_cast<IntType>(1) << bitPosition;
    }
}

template <typename Derived, size_t N>
void BitArrayBase<Derived, N>::set(size_t index, bool value)
{
    if (value)
        set(index);
    else
        reset(index);
}

// Reset Bit Method Implementation
template <typename Derived, size_t N>
void BitArrayBase<Derived, N>::reset(size_t index)
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
        size_t arrayIndex = index / szInt;
        size_t bitPosition = index % szInt;
        bits_[arrayIndex] &= ~(static_cast<IntType>(1) << bitPosition);
    }
}

// Right Shift Operator Implementation
template <typename Derived, size_t N>
Derived &BitArrayBase<Derived, N>::operator>>=(size_t shift)
{
    if (shift >= N)
    {
        std::fill(std::begin(bits_), std::end(bits_), 0); // Clear all bits if shift >= N
        return *static_cast<Derived *>(this);
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
    return *static_cast<Derived *>(this);
}

// Left Shift Operator Implementation
template <typename Derived, size_t N>
Derived &BitArrayBase<Derived, N>::operator<<=(size_t shift)
{
    if (shift >= N)
    {
        std::fill(std::begin(bits_), std::end(bits_), 0); // Clear all bits if shift >= N
        return *static_cast<Derived *>(this);
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
    return *static_cast<Derived *>(this);
}

// Bitwise OR Assignment Operator Implementation
template <typename Derived, size_t N>
Derived &BitArrayBase<Derived, N>::operator|=(const Derived &other)
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
    return *static_cast<Derived *>(this);
}

// Bitwise AND Assignment Operator Implementation
template <typename Derived, size_t N>
Derived &BitArrayBase<Derived, N>::operator&=(const Derived &other)
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
    return *static_cast<Derived *>(this);
}

// Bitwise XOR Assignment Operator Implementation
template <typename Derived, size_t N>
Derived &BitArrayBase<Derived, N>::operator^=(const Derived &other)
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
    return *static_cast<Derived *>(this);
}

// Unary NOT Operator Implementation
template <typename Derived, size_t N>
Derived BitArrayBase<Derived, N>::operator~() const
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
    constexpr size_t excessBits = N % szInt;
    if (excessBits != 0)
    {
        result.bits_[numInt - 1] &= (static_cast<IntType>(1) << excessBits) - 1;
    }
    return result;
}

// Comparison Operator Implementations
template <typename Derived, size_t N>
bool BitArrayBase<Derived, N>::operator<(const Derived &other) const
{
    if constexpr (numInt == 1)
    {
        return (bits_[0] < other.bits_[0]);
    }
    else
    {
        for (size_t i = numInt; i > 0; --i)
        {
            if (bits_[i - 1] < other.bits_[i - 1])
                return true;
        }
    }
    return false;
}

template <typename Derived, size_t N>
bool BitArrayBase<Derived, N>::operator>=(const Derived &other) const
{
    return !(*this < other);
}

template <typename Derived, size_t N>
bool BitArrayBase<Derived, N>::operator>(const Derived &other) const
{
    return other < *this;
}

template <typename Derived, size_t N>
bool BitArrayBase<Derived, N>::operator<=(const Derived &other) const
{
    return !(other < *this);
}

template <typename Derived, size_t N>
bool BitArrayBase<Derived, N>::operator==(const Derived &other) const
{
    if constexpr (numInt == 1)
    {
        return (bits_[0] == other.bits_[0]);
    }
    else
    {
        for (size_t i = 0; i < numInt; ++i)
        {
            if (bits_[i] != other.bits_[i])
                return false;
        }
    }
    return true;
}

template <typename Derived, size_t N>
bool BitArrayBase<Derived, N>::operator!=(const Derived &other) const
{
    return !(*this == other);
}

// All bits check implementation
template <typename Derived, size_t N>
bool BitArrayBase<Derived, N>::all()
{
    if constexpr (numInt == 1)
    {
        if constexpr (N % szInt != 0)
        {
            return ((bits_[0] | ((static_cast<IntType>(-1) << (N % szInt)))) == static_cast<IntType>(-1));
        }
        else
        {
            return bits_[0] == static_cast<IntType>(-1);
        }
    }
    else
    {
        for (auto i(0); i < numInt - 1; ++i)
        {
            if (bits_[i] != static_cast<IntType>(-1))
                return false;
        }
        if constexpr (N % szInt != 0)
        {
            return ((bits_[numInt - 1] | ((static_cast<IntType>(-1) << (N % szInt)))) == static_cast<IntType>(-1));
        }
        else
        {
            return bits_[numInt - 1] == static_cast<IntType>(-1);
        }
    }
}

// Any bit set check implementation
template <typename Derived, size_t N>
bool BitArrayBase<Derived, N>::any()
{
    if constexpr (numInt == 1)
    {
        return (bits_[0] != 0);
    }
    else
    {
        for (auto i(0); i < numInt; ++i)
        {
            if (bits_[i] != 0)
                return true;
        }
        return false;
    }
}

// No bits set check implementation
template <typename Derived, size_t N>
bool BitArrayBase<Derived, N>::none()
{
    return !any();
}

// Counting consecutive zeros form the right implementation
template <typename Derived, size_t N>
size_t BitArrayBase<Derived, N>::countr_zero() const
{
    size_t cnt = std::countr_zero(bits_[0]);
    if constexpr (numInt == 1)
    {       
        if constexpr (N < szInt)
        {
            cnt = (cnt>N ? N : cnt);
        }
    }
    else
    {
        if (cnt != szInt){
             return cnt;
        } 
        else
        {
            int cntInt;
            for (size_t i = 1; i < numInt; ++i)
            {
                cntInt = std::countr_zero(bits_[i]);
                cnt += cntInt;
                if (cntInt != szInt){
                    return cnt;
                } 
            }
            if constexpr (N % szInt != 0)
            {
                cnt = (cnt>N ? N : cnt);
            }
        }
    }
    return cnt;
}

// BitReference support implementation
template <typename Derived, size_t N>
BitReference<Derived, N> BitArrayBase<Derived, N>::operator[](size_t index)
{
    return BitReference<Derived, N>(*this, index);
}

// Stream output operator implementation
template <typename Derived, size_t N>
std::ostream &operator<<(std::ostream &os, const BitArrayBase<Derived, N> &bitArray)
{
    for (size_t i = N; i > 0; --i)
    {
        os << bitArray.test(i - 1); // Output each bit
    }
    return os;
}

// Right shift operator
template <typename Derived, size_t N>
Derived operator>>(const BitArrayBase<Derived, N> &bitArray, size_t shift) {
    Derived result(static_cast<const Derived&>(bitArray)); // Cast to Derived type and copy construct
    result >>= shift; // Apply right shift using member operator>>=
    return result; // Return the result
}

// Left shift operator
template <typename Derived, size_t N>
Derived operator<<(const BitArrayBase<Derived, N> &bitArray, size_t shift) {
    Derived result(static_cast<const Derived&>(bitArray)); // Cast to Derived type and copy construct
    result <<= shift; // Apply left shift using member operator<<=
    return result; // Return the result
}

// Bitwise OR operator
template <typename Derived, size_t N>
Derived operator|(const BitArrayBase<Derived, N> &bitArray1, const BitArrayBase<Derived, N> &bitArray2) {
    Derived result(static_cast<const Derived&>(bitArray1)); // Cast to Derived type and copy construct
    result |= static_cast<const Derived&>(bitArray2); // Apply bitwise OR using member operator|=
    return result; // Return the result
}

// Bitwise AND operator
template <typename Derived, size_t N>
Derived operator&(const BitArrayBase<Derived, N> &bitArray1, const BitArrayBase<Derived, N> &bitArray2) {
    Derived result(static_cast<const Derived&>(bitArray1)); // Cast to Derived type and copy construct
    result &= static_cast<const Derived&>(bitArray2); // Apply bitwise AND using member operator&=
    return result; // Return the result
}

// Bitwise XOR operator
template <typename Derived, size_t N>
Derived operator^(const BitArrayBase<Derived, N> &bitArray1, const BitArrayBase<Derived, N> &bitArray2) {
    Derived result(static_cast<const Derived&>(bitArray1)); // Cast to Derived type and copy construct
    result ^= static_cast<const Derived&>(bitArray2); // Apply bitwise XOR using member operator^=
    return result; // Return the result
}

//Bitreference contructor
template <typename Derived, size_t N>
BitReference<Derived, N>::BitReference(BitArrayBase<Derived, N> &ba, size_t idx): bitArray(ba), index(idx) {}

// Assignment operator to set the bit
template <typename Derived, size_t N>
BitReference<Derived, N> &BitReference<Derived, N>::operator=(bool value)
{
    bitArray.set(index, value);
    return *this;
}

// Conversion to bool (read the bit)
template <typename Derived, size_t N>
BitReference<Derived, N>::operator bool() const
{
    return bitArray.test(index);
}


#endif // BITARRAY_HPP
