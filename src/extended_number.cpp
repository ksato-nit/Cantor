#include "extended_number.hpp"

ExtendedNumber::ExtendedNumber(){
    this->re = Number::ZERO();
    this->im = Number::ZERO();
}

ExtendedNumber::ExtendedNumber(int c){
    this->re = Number(c);
    this->im = Number::ZERO();
}

ExtendedNumber::ExtendedNumber(Number re, Number im){
    this->re = re;
    this->im = im;
}

ExtendedNumber ExtendedNumber::operator + (const ExtendedNumber& c) const{
    ExtendedNumber z(this->re + c.re, this->im + c.im);
    return z;
}

ExtendedNumber ExtendedNumber::operator - (const ExtendedNumber& c) const{
    ExtendedNumber z(this->re - c.re, this->im - c.im);
    return z;
}

ExtendedNumber ExtendedNumber::operator * (const ExtendedNumber& c) const{
    Number x = this->re * c.re - this->im * c.im;
    Number y = this->re * c.im + this->im * c.re;
    ExtendedNumber z(x, y);
    return z;
}

ExtendedNumber ExtendedNumber::operator * (const Number m) const{
    ExtendedNumber z(this->re * m, this->im * m);
    return z;
}

ExtendedNumber ExtendedNumber::operator * (int m) const{
    ExtendedNumber z(this->re * m, this->im * m);
    return z;
}

ExtendedNumber ExtendedNumber::operator / (const ExtendedNumber& c) const{
    ExtendedNumber cinv = c.inv();
    return (*this) * cinv;
}

bool ExtendedNumber::operator == (const ExtendedNumber& c) const{
    return (this->re - c.re).isZero() && (this->im - c.im).isZero();
}

bool ExtendedNumber::operator != (const ExtendedNumber& c) const{
    return !(*this == c);
}

ExtendedNumber ExtendedNumber::operator + () const{
    return *this;
}

ExtendedNumber ExtendedNumber::operator - () const{
    ExtendedNumber m(this->re * Number::MINUS_ONE(), this->im * Number::MINUS_ONE());
    return m;
}

ExtendedNumber ExtendedNumber::ZERO(){
    ExtendedNumber zero(0, 0);
    return zero;
}

ExtendedNumber ExtendedNumber::ONE(){
    ExtendedNumber one(1, 0);
    return one;
}

ExtendedNumber ExtendedNumber::MINUS_ONE(){
    ExtendedNumber minus_one(-1, 0);
    return minus_one;
}

ExtendedNumber ExtendedNumber::inv() const{
    Number denom = this->re * this->re + this->im * this->im;
    denom = denom.inv();
    Number x = this->re * denom;
    Number y = -this->im * denom;
    ExtendedNumber z(x, y);
    return z;
}

bool ExtendedNumber::isZero() const{
    return this->re.isZero() && this->im.isZero();;
}
