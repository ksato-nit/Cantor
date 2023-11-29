#include "extended_number.hpp"

ExtendedNumber::ExtendedNumber(Polynomial f){
    this->g.resize(1);
    this->g.coeff[0] = Number::ZERO();
    this->g.coeff[1] = Number::ZERO();
    this->f = f;
}

ExtendedNumber::ExtendedNumber(Polynomial f, Number x, Number y){
    this->g.resize(1);
    this->g.coeff[0] = x;
    this->g.coeff[1] = y;
    this->f = f;
}

ExtendedNumber::ExtendedNumber(Polynomial f, Polynomial g){
    this->g = g % f;
    this->f = f;
}

ExtendedNumber ExtendedNumber::operator + (const ExtendedNumber& c) const{
    Polynomial h = (this->g + c.g) % this->f;
    ExtendedNumber z(f, h);
    return z;
}

ExtendedNumber ExtendedNumber::operator - (const ExtendedNumber& c) const{
    Polynomial h = (this->g * c.g) % this->f;
    ExtendedNumber z(f, h);
    return z;
}

ExtendedNumber ExtendedNumber::operator * (const ExtendedNumber& c) const{
    Polynomial h = (this->g * c.g) % this->f;
    ExtendedNumber z(f, h);
    return z;
}

ExtendedNumber ExtendedNumber::operator * (const Number m) const{
    Polynomial h = (this->g * m) % this->f;
    ExtendedNumber z(f, h);
    return z;
}

ExtendedNumber ExtendedNumber::operator / (const ExtendedNumber& c) const{
    ExtendedNumber cinv = c.inv();
    return (*this) * cinv;
}

bool ExtendedNumber::operator == (const ExtendedNumber& c) const{
    return (this->g - c.g).isZero();
}

bool ExtendedNumber::operator != (const ExtendedNumber& c) const{
    return !(*this == c);
}

ExtendedNumber ExtendedNumber::operator + () const{
    return *this;
}

ExtendedNumber ExtendedNumber::operator - () const{
    ExtendedNumber m(this->f, this->g * Number::MINUS_ONE());
    return m;
}

ExtendedNumber ExtendedNumber::ZERO(){
    ExtendedNumber zero(1, 0);
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
    auto tup = Polynomial::extended_gcd(this->g, this->f);
    Polynomial h = std::get<1>(tup);
    ExtendedNumber z(this->f, h);
    return z;
}

bool ExtendedNumber::isZero(){
    return this->g.isZero();
}
