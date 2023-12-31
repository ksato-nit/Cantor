#pragma once
#include "algorithm"
#include "iostream"
#include "string"
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/integer/mod_inverse.hpp>

// 抽象的な体の要素を実現する．
class Number {
public:
    static boost::multiprecision::int1024_t CHARA;

    boost::multiprecision::int1024_t value;

    Number();
    Number(int);
    Number(boost::multiprecision::int1024_t);
    Number(const Number&);
    Number(Number&&) noexcept;
    Number operator + (const Number&) const;
    Number operator - (const Number&) const;
    Number operator * (const Number&) const;
    Number operator * (const int) const;
    Number operator * (const boost::multiprecision::int1024_t) const;
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
