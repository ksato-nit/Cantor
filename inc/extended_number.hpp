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
    Polynomial f; // 既約多項式

    Polynomial g;

    ExtendedNumber(Polynomial);
    ExtendedNumber(Polynomial, Polynomial);
    ExtendedNumber(Polynomial, Number, Number);

    ExtendedNumber operator + (const ExtendedNumber&) const;
    ExtendedNumber operator - (const ExtendedNumber&) const;
    ExtendedNumber operator * (const ExtendedNumber&) const;
    ExtendedNumber operator * (const Number) const;
    ExtendedNumber operator / (const ExtendedNumber&) const;
    ExtendedNumber inv() const;

    bool operator == (const ExtendedNumber& c) const;
    bool operator != (const ExtendedNumber& c) const;
    ExtendedNumber operator + () const;
    ExtendedNumber operator - () const;

    static ExtendedNumber ZERO();
    static ExtendedNumber ONE();
    static ExtendedNumber MINUS_ONE();

    bool isZero();
};

inline std::ostream& operator << (std::ostream& os, const ExtendedNumber& n){
    os << n.g;
    return os;
}
