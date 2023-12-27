#include "number.hpp"

NTL::ZZ Number::CHARA;

Number::Number(){
}

Number::Number(NTL::ZZ_p x){
    this->value = x;
}

Number::Number(int x){
    this->value = x;
}

Number::Number(const Number& num){
    this->value = num.value;
}

Number::Number(Number&& y) noexcept{
    this->value = std::move(y.value);
}

void Number::set_str(const char* str, const int base){
    this->value = NTL::conv<NTL::ZZ_p>(str);
    return;
}

Number Number::operator + (const Number& y) const{
    Number z;
    z.value = this->value + y.value;
    return z;
}

Number Number::operator - (const Number& y) const{
    Number z;
    z.value = this->value - y.value;
    return z;
}

Number Number::operator * (const Number& y) const{
    Number z;
    z.value = (this->value) * (y.value);
    return z;
}

Number Number::operator * (const int y) const{
    Number z;
    z.value = this->value * y;
    return z;
}

Number Number::operator * (const NTL::ZZ_p y) const{
    Number z;
    z.value = this->value * y;
    return z;
}

Number Number::operator / (const Number& y) const{
    return (*this) * y.inv();
}

bool Number::operator == (const Number& y) const{
    NTL::ZZ_p ev;
    ev = this->value - y.value;
    return (ev == 0);
}

bool Number::operator != (const Number& y) const{
    return !(*this == y);
}

Number Number::operator + () const{
    return *this;
}

Number Number::operator - () const{
    Number z;
    z.value = -this->value;
    return z;
}

Number& Number::operator = (const Number& y) {
    this->value = y.value;
    return *this;
}

Number& Number::operator = (Number&& y) noexcept{
    this->value = std::move(y.value);
    return *this;
}

Number Number::ZERO(){
    Number ret(0);
    return ret;
}

Number Number::ONE(){
    Number ret(1);
    return ret;
}

Number Number::MINUS_ONE(){
    Number ret(-1);
    return ret;
}

Number Number::inv() const{
    Number ret;
    ret.value = NTL::inv(this->value);
    return ret;
}

bool Number::isZero() const{
    return (this->value == 0);
}
