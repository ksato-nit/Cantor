#pragma once
#include "algorithm"
#include "iostream"
#include "string"
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>

// 抽象的な体の要素を実現する．
class Number {
public:
    static NTL::ZZ CHARA;

    NTL::ZZ_p value;

    Number();
    Number(int);
    Number(NTL::ZZ_p);
    Number(const Number&);
    Number(Number&&) noexcept;
    Number operator + (const Number&) const;
    Number operator - (const Number&) const;
    Number operator * (const Number&) const;
    Number operator * (const int) const;
    Number operator * (const NTL::ZZ_p) const;
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
