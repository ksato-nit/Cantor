#include "number.hpp"

const mpz_class Number::CHARA("16613960161207197506610974848157905611744466278275346794947826509160636299163", 10);

Number::Number(){
    this->value = 0;
}

Number::Number(mpz_class x){
    this->value = x;
}

Number Number::operator + (const Number& y) const{
    Number z((y.value + this->value) % CHARA);
    return z;
}

Number Number::operator - (const Number& y) const{
    Number z((this->value - y.value) % CHARA);
    return z;
}

Number Number::operator * (const Number& y) const{
    Number z((y.value * this->value) % CHARA);
    return z;
}

Number Number::operator * (const mpz_class y) const{
    Number z((y * this->value) % CHARA);
    return z;
}

Number Number::operator / (const Number& y) const{
    return (*this) * y.inv();
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
    mpz_class value = -this->value;
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

    mpz_class s = this->value, t = CHARA;
    mpz_class x = 1, u = 0;
    mpz_class k;
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
