#include <bitset>
#include <iostream>

using std::bitset;

template<size_t Dim, size_t N>
class Morton: public bitset<Dim*N> {
public:
    using bitset<Dim*N>::bitset;

    Morton(const bitset<Dim*N>& bs) : bitset<Dim*N>(bs) {}

    Morton<Dim, N>& operator<<(size_t pos)
    {
        bitset<Dim*N>& base = *this;
        base = base << (Dim*pos);
        return *this;
    }

    Morton<Dim, N>& operator>>(size_t pos)
    {
        bitset<Dim*N>& base = *this;
        base = base >> (Dim*pos);
        return *this;
    }

    // Addition operators
    Morton operator+(const Morton &other) const
    {
        auto carry = (*this) & other;  // Initial carry
        auto result = (*this) ^ other; // Initial sum (without carry)
        // Propagate carry
        while (carry.any())
        {
            auto shiftedCarry = carry << 1;
            carry = result & shiftedCarry;
            result ^= shiftedCarry;
        }
        return result;
    }

    Morton& operator+=(const Morton& other) {
        auto carry = (*this) & other;  // Initial carry
        (*this) ^= other; // Initial sum (without carry)
        // Propagate carry
        while (carry.any())
        {
            auto shiftedCarry = carry << 1;
            carry = (*this) & shiftedCarry;
            (*this) ^= shiftedCarry;
        }
        return *this;
    }

    // Subtraction operators
    Morton operator-(const Morton& other) const {
        auto complement = ~other;
        complement += ones;
        return (*this) + complement;
    }

    Morton& operator-=(const Morton& other) {
        auto complement = ~other;
        complement += ones;
        (*this) += complement;
        return *this;
    }
/*
    bool operator<(const Morton& other) const {
       return value < other.value;
    }
    bool operator>(const Morton& other) const {
        return value > other.value;
    }
    bool operator<=(const Morton& other) const {
        return value <= other.value;
    }
    bool operator>=(const Morton& other) const {
        return value >= other.value;
    }
*/
    private:

    static constexpr Morton<Dim, N> createOnes() {
        Morton<Dim, N> result;
        for (size_t i = 0; i < Dim; ++i) {
            result.set(i);
        }
        return result;
    }
    public:
    static constexpr Morton<Dim,N> ones = createOnes();
};

int main() {
    Morton<2,21> num1(1);
    Morton<2,21> num2(0);

    num1 << 1;
    num2 +=1;
    num1 -=2;
    std::cout << "num1: " << num1 << ", num2: " << num2 << std::endl;

    return 0;
}
