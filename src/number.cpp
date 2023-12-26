#include "number.hpp"

mpz_class Number::CHARA;
mpz_class Number::MCHARA;

const mpz_class Number::zero = 0;
const mpz_class Number::one = 1;
const mpz_class Number::minus_one = -1;

Number::Number(){
}

Number::Number(mpz_class x){
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
    this->value.set_str(str, base);
    return;
}

Number Number::operator + (const Number& y) const{
    Number z;
    z.value = this->value + y.value;
    if(z.value > CHARA){
        z.value -= CHARA;
    }else if(z.value < MCHARA){
        z.value += CHARA;
    }
    return z;
}

Number Number::operator - (const Number& y) const{
    Number z;
    z.value = this->value - y.value;
    if(z.value > CHARA){
        z.value -= CHARA;
    }else if(z.value < MCHARA){
        z.value += CHARA;
    }
    return z;
}

Number Number::operator * (const Number& y) const{
    Number z;
    z.value = this->value * y.value;
    z.value %= CHARA;
    return z;
}

Number Number::operator * (const int y) const{
    Number z;
    z.value = this->value * y;
    return z;
}

Number Number::operator * (const mpz_class y) const{
    Number z;
    z.value = this->value * y;
    return z;
}

Number Number::operator / (const Number& y) const{
    return (*this) * y.inv();
}

bool Number::operator == (const Number& y) const{
    mpz_class ev;
    ev = this->value - y.value;
    ev = ev % CHARA;
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
    mpz_t inv;
    mpz_init(inv);
    mpz_invert(inv, this->value.get_mpz_t(), CHARA.get_mpz_t());
    Number ret;
    ret.value = mpz_class(inv);
    return ret;
}

bool Number::isZero() const{
    return (this->value % CHARA == 0);
}
