#pragma once
#include "algorithm"
#include "iostream"

// 抽象的な体の要素を実現する．
class Number {
public:
    static const int CHARA = 31;
        
    int value;

    Number();
    Number(int value);
    Number operator + (const Number&) const;
    Number operator - (const Number&) const;
    Number operator * (const Number&) const;
    Number operator * (const int ) const;
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
