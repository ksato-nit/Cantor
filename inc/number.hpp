#pragma once
#include "algorithm"
#include "iostream"
#include <boost/multiprecision/cpp_int.hpp>
#include <gmpxx.h>

// 抽象的な体の要素を実現する．
class Number {
public:
    mpz_class CHARA;
    mpz_class value;

    Number();
    Number(mpz_class value);
    Number operator + (const Number&) const;
    Number operator - (const Number&) const;
    Number operator * (const Number&) const;
    Number operator * (const mpz_class) const;
    Number operator / (const Number&) const;
    bool operator == (const Number&) const;
    bool operator != (const Number&) const;
    Number operator + () const;
    Number operator - () const;
    Number inv() const;
    bool isZero();

    static Number ZERO();
    static Number ONE();
    static Number MINUS_ONE();
};

inline std::ostream& operator << (std::ostream& os, const Number& x){
    os << x.value;
    return os;
}
