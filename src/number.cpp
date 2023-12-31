#include "number.hpp"

boost::multiprecision::int1024_t Number::CHARA;

Number::Number(){
}

Number::Number(boost::multiprecision::int1024_t x){
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
    this->value.assign(str);
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

Number Number::operator * (const boost::multiprecision::int1024_t y) const{
    Number z;
    z.value = this->value * y;
    return z;
}

Number Number::operator / (const Number& y) const{
    return (*this) * y.inv();
}

bool Number::operator == (const Number& y) const{
    boost::multiprecision::int1024_t ev;
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
    if(this->value < 0){
        Number mx(this->value * -1);
        mx = mx.inv();
        mx.value *= -1;
        return mx;
    }

    boost::multiprecision::int1024_t s = this->value, t = CHARA;
    boost::multiprecision::int1024_t x = 1, u = 0;
    boost::multiprecision::int1024_t k;
    while(t > 0) {
        k = s / t;
        s = s - (k * t);
        std::swap(s, t);
        x = x - (k * u);
        std::swap(x, u);
    }
    x = x % CHARA;
    if(x < 0){
        x = x + CHARA;
    }
    return x;
}

bool Number::isZero() const{
    return (this->value == 0);
}
