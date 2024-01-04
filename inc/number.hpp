#pragma once
#include "algorithm"
#include "iostream"
#include "string"
#include <boost/multiprecision/cpp_int.hpp>
#include <gmpxx.h>

// 抽象的な体の要素を実現する．
class Number {
public:
    static mpz_t CHARA;
    static mpz_t MCHARA; // -CHARA

    mpz_t value;

    Number();
    ~Number();
    Number(int);
    Number(const Number&);
    Number(Number&&) noexcept;
    Number operator + (const Number&) const;
    Number operator - (const Number&) const;
    Number operator * (const Number&) const;
    Number operator * (const int) const;
    Number operator * (const mpz_t) const;
    Number operator / (const Number&) const;
    void operator += (const Number&);
    void operator -= (const Number&);
    void operator *= (const Number&);
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
