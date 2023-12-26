#pragma once
#include "algorithm"
#include "iostream"
#include "string"
#include <boost/multiprecision/cpp_int.hpp>
#include <gmpxx.h>

// 抽象的な体の要素を実現する．
class Number {
public:
    static mpz_class CHARA;
    static mpz_class MCHARA; // -CHARA

    static const mpz_class zero;
    static const mpz_class one;
    static const mpz_class minus_one;

    mpz_class value;

    Number();
    Number(int);
    Number(mpz_class);
    Number(const Number&);
    Number(Number&&) noexcept;
    Number operator + (const Number&) const;
    Number operator - (const Number&) const;
    Number operator * (const Number&) const;
    Number operator * (const int) const;
    Number operator * (const mpz_class) const;
    Number operator / (const Number&) const;
    bool operator == (const Number&) const;
    bool operator != (const Number&) const;
    Number operator + () const;
    Number operator - () const;
    Number& operator = (const Number &);
    Number& operator = (Number &&) noexcept;
    Number inv() const;
    bool isZero() const;
    void set_str(const char*, const int);

    static Number ZERO();
    static Number ONE();
    static Number MINUS_ONE();
};

inline std::ostream& operator << (std::ostream& os, const Number& x){
    os << x.value;
    return os;
}
