#pragma once
#include "algorithm"
#include "iostream"
#include "number.hpp"
#include <gmpxx.h>

// 2 次拡大体の要素．
class ExtendedNumber{
public:
    static mpz_t CHARA;
    mpz_t re, im;
    ExtendedNumber();
    ExtendedNumber(int, int);
    ExtendedNumber(mpz_t, mpz_t);

    ExtendedNumber operator * (const ExtendedNumber&) const;

    bool isZero();
    
    static ExtendedNumber ZERO();
    static ExtendedNumber ONE();
};

inline std::ostream& operator << (std::ostream& os, const ExtendedNumber& n){
    os << n.re << " + i" << n.im;
    return os;
}
