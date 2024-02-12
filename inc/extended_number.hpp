#pragma once
#include "algorithm"
#include "iostream"
#include "number.hpp"
#include "polynomial.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <gmpxx.h>

// 2 次拡大体の要素．
class ExtendedNumber{
public:
    Number re, im;

    ExtendedNumber();
    ExtendedNumber(int);
    ExtendedNumber(Number, Number);

    ExtendedNumber operator + (const ExtendedNumber&) const;
    ExtendedNumber operator - (const ExtendedNumber&) const;
    ExtendedNumber operator * (const ExtendedNumber&) const;
    ExtendedNumber operator * (const Number) const;
    ExtendedNumber operator * (const int) const;
    ExtendedNumber operator / (const ExtendedNumber&) const;
    ExtendedNumber inv() const;

    bool operator == (const ExtendedNumber& c) const;
    bool operator != (const ExtendedNumber& c) const;
    ExtendedNumber operator + () const;
    ExtendedNumber operator - () const;

    static ExtendedNumber ZERO();
    static ExtendedNumber ONE();
    static ExtendedNumber MINUS_ONE();

    bool isZero() const;
};

inline std::ostream& operator << (std::ostream& os, const ExtendedNumber& n){
    os << n.re << " + i" << n.im;
    return os;
}
