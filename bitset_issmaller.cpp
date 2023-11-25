#include <bitset>
#include <iostream>
#include <type_traits>
#include <limits>

template <size_t N>
struct IsSmaller
{
    bool operator()(const std::bitset<N> &a, const std::bitset<N> &b) const
    {
        if constexpr (N <= sizeof(unsigned long) * 8)
        {
            // If the bitset fits into an unsigned long, use to_ulong for comparison
            return a.to_ulong() < b.to_ulong();
        }
        else if constexpr (N <= sizeof(unsigned long long) * 8)
        {
            // If the bitset fits into an unsigned long long, use to_ullong for comparison
            return a.to_ullong() < b.to_ullong();
        }
        else
        {
            constexpr unsigned long max_ulong = std::numeric_limits<unsigned long>::max();
            constexpr std::bitset<N> mask(max_ulong);
            constexpr size_t ulong_bits = sizeof(unsigned long) * 8;
            constexpr size_t iterations = (N + ulong_bits - 1) / ulong_bits; // Calculate number of segments
            unsigned long segments_a[iterations], segments_b[iterations];
            std::bitset<N> temp_a = a, temp_b = b;
            for (size_t i = 0; i < iterations; ++i)
            {
                segments_a[i] = (temp_a & mask).to_ulong();
                temp_a >>= ulong_bits;
                segments_b[i] = (temp_b & mask).to_ulong();
                temp_b >>= ulong_bits;
            }
            for (size_t i = iterations - 1; i != (-1); --i)
                if (segments_a[i] != segments_b[i])
                {
                    return segments_a[i] < segments_b[i];
                }
            return false; // Equal if all segments are the same
        }
    }
};

int main()
{
    std::bitset<128> bitset1(1);
    std::bitset<128> bitset2(1);
    IsSmaller<128> isSmaller;

    bitset1 <<= 111;
    bitset1 |= 1;
    bitset2 <<= 111;
    bitset2 |= 2;

    if (isSmaller(bitset1, bitset2))
    {
        std::cout << "bitset1 is smaller" << std::endl;
    }
    else
    {
        std::cout << "bitset2 is smaller or they are equal" << std::endl;
    }

    return 0;
}