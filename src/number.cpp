#include "number.hpp"

mpz_t Number::CHARA;
mpz_t Number::MCHARA;

Number::Number(){
    mpz_init2(this->value, 512);
}

Number::~Number(){
    mpz_clear(this->value);
}

Number::Number(int x){
    mpz_init_set_si(this->value, x);
}

Number::Number(const Number& num){
    mpz_init2(this->value, 512);
    mpz_set(this->value, num.value);
}

Number::Number(Number&& y) noexcept{
    this->value->_mp_alloc = y.value->_mp_alloc;
    this->value->_mp_size = y.value->_mp_size;
    this->value->_mp_d = y.value->_mp_d;
    mpz_init_set_ui(y.value, 0);
}

void Number::set_str(const char* str, const int base){
    mpz_set_str(this->value, str, base);
    return;
}

void Number::operator += (const Number& y){
    mpz_add(this->value, this->value, y.value);
    if(mpz_cmp(this->value, CHARA) >= 0){
        mpz_sub(this->value, this->value, CHARA);
    }else if(mpz_cmp(this->value, MCHARA) <= 0){
        mpz_add(this->value, this->value, CHARA);
    }
    return;
}

void Number::operator -= (const Number& y){
    mpz_sub(this->value, this->value, y.value);
    if(mpz_cmp(this->value, CHARA) >= 0){
        mpz_sub(this->value, this->value, CHARA);
    }else if(mpz_cmp(this->value, MCHARA) <= 0){
        mpz_add(this->value, this->value, CHARA);
    }
    return;
}

void Number::operator *= (const Number& y){
    mpz_mul(this->value, this->value, y.value);
    mpz_mod(this->value, this->value, CHARA);
    return;
}

Number Number::operator + (const Number& y) const{
    Number z;
    mpz_add(z.value, this->value, y.value);
    //mpz_mod(z.value, z.value, CHARA);
    if(mpz_cmp(z.value, CHARA) >= 0){
        mpz_sub(z.value, z.value, CHARA);
    }else if(mpz_cmp(z.value, MCHARA) <= 0){
        mpz_add(z.value, z.value, CHARA);
    }
    return z;
}

Number Number::operator - (const Number& y) const{
    Number z;
    mpz_sub(z.value, this->value, y.value);
    //mpz_mod(z.value, z.value, CHARA);
    if(mpz_cmp(z.value, CHARA) >= 0){
        mpz_sub(z.value, z.value, CHARA);
    }else if(mpz_cmp(z.value, MCHARA) <= 0){
        mpz_add(z.value, z.value, CHARA);
    }
    return z;
}

Number Number::operator * (const Number& y) const{
    Number z;
    mpz_mul(z.value, this->value, y.value);
    mpz_mod(z.value, z.value, CHARA);
    return z;
}

Number Number::operator * (const int y) const{
    Number z;
    mpz_t ym;
    mpz_init_set_ui(ym, y);
    mpz_mul(z.value, this->value, ym);
    mpz_mod(z.value, z.value, CHARA);
    return z;
}

Number Number::operator * (const mpz_t y) const{
    Number z;
    mpz_mul(z.value, this->value, y);
    mpz_mod(z.value, z.value, CHARA);
    return z;
}

Number Number::operator / (const Number& y) const{
    return (*this) * y.inv();
}

bool Number::operator == (const Number& y) const{
    mpz_t ev;
    mpz_init(ev);
    mpz_sub(ev, this->value, y.value);
    mpz_mod(ev, ev, CHARA);
    return (mpz_sgn(ev) == 0);
}

bool Number::operator != (const Number& y) const{
    return !(*this == y);
}

Number Number::operator + () const{
    return *this;
}

Number Number::operator - () const{
    Number z;
    mpz_neg(z.value, this->value);
    return z;
}

Number& Number::operator = (const Number& y) {
    mpz_set(this->value, y.value);
    return *this;
}

Number& Number::operator = (Number&& y) noexcept{
    this->value->_mp_alloc = y.value->_mp_alloc;
    this->value->_mp_size = y.value->_mp_size;
    this->value->_mp_d = y.value->_mp_d;
    mpz_init_set_ui(y.value, 0);
    return *this;
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
    Number ret;
    mpz_invert(ret.value, this->value, CHARA);

    mpz_mod(ret.value, ret.value, CHARA);

    return ret;
}

bool Number::isZero() const{
    return (mpz_sgn(this->value) == 0);
}
