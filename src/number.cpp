#include "number.hpp"

Number::Number(){
    this->value = 0;
}

Number::Number(int x){
    this->value = x % CHARA;
}

Number Number::operator + (const Number& y) const{
    int value = (y.value + this->value) % CHARA;
    Number z(value);
    return z;
}

Number Number::operator - (const Number& y) const{
    int value = (this->value - y.value) % CHARA;
    Number z(value);
    return z;
}

Number Number::operator * (const Number& y) const{
    int value = (y.value * this->value) % CHARA;
    Number z(value);
    return z;
}

Number Number::operator * (const int y) const{
    int value = (y * this->value) % CHARA;
    Number z(value);
    return z;
}

Number Number::operator / (const Number& y) const{
    Number yinv = y.inv();
    int value = (this->value * yinv.value) % CHARA;
    Number z(value);
    return z;
}

bool Number::operator == (const Number& y) const{
    return ((this->value - y.value) % CHARA) == 0;
}

bool Number::operator != (const Number& y) const{
    return !(*this == y);
}

Number Number::operator + () const{
    return *this;
}

Number Number::operator - () const{
    int value = -this->value % CHARA;
    Number z(value);
    return z;
}

Number Number::ZERO(){
    Number zero(0);
    return zero;
}

Number Number::ONE(){
    Number one(1);
    return one;
}

Number Number::MINUS_ONE(){
    Number minus_one(-1);
    return minus_one;
}

Number Number::inv() const{
    if(this->value < 0){
        Number mx(this->value * -1);
        mx = mx.inv();
        mx.value *= -1;
        return mx;
    }

    int s = this->value, t = CHARA;
    int x = 1, u = 0;
    int k;
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

bool Number::isZero(){
    return this->value == 0;
}
