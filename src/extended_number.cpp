#include "extended_number.hpp"

mpz_t ExtendedNumber::CHARA;
mpz_t ExtendedNumber::MCHARA;

ExtendedNumber::ExtendedNumber(){
    mpz_init_set_ui(this->re, 0);
    mpz_init_set_ui(this->im, 0);
}

ExtendedNumber::ExtendedNumber(int re, int im){
    mpz_init_set_ui(this->re, re);
    mpz_init_set_ui(this->im, im);
}

ExtendedNumber::ExtendedNumber(mpz_t re, mpz_t im){
    mpz_init_set(this->re, re);
    mpz_init_set(this->im, im);
}

ExtendedNumber ExtendedNumber::operator * (const ExtendedNumber& y) const{
    ExtendedNumber z;
    mpz_t temp;
    mpz_init(temp);
    mpz_mul(z.re, this->re, y.re);
    mpz_mod(z.re, z.re, CHARA);
    mpz_mul(temp, this->im, y.im);
    mpz_mod(temp, temp, CHARA);
    mpz_sub(z.re, z.re, temp);
    mpz_mul(z.im, this->re, y.im);
    mpz_mod(z.im, z.im, CHARA);
    mpz_mul(temp, this->im, y.re);
    mpz_mod(temp, temp, CHARA);
    mpz_sub(z.im, z.im, temp);
    return z;
}

bool ExtendedNumber::operator == (const ExtendedNumber& y) const{
    mpz_t ev;
    mpz_init(ev);
    mpz_sub(ev, this->re, y.re);
    mpz_mod(ev, ev, CHARA);
    if(mpz_sgn(ev) != 0) return false;

    mpz_sub(ev, this->im, y.im);
    mpz_mod(ev, ev, CHARA);
    return (mpz_sgn(ev) == 0);
}

bool ExtendedNumber::operator != (const ExtendedNumber& y) const{
    return !(*this == y);
}

bool ExtendedNumber::isZero() const{
    mpz_t ev;
    mpz_init_set(ev, this->re);
    mpz_mod(ev, ev, CHARA);
    if(mpz_sgn(ev) != 0) return false;
    mpz_set(ev, this->im);
    mpz_mod(ev, ev, CHARA);
    return (mpz_sgn(ev) == 0);
}

ExtendedNumber ExtendedNumber::ZERO(){
    ExtendedNumber zero(0, 0);
    return zero;
}

ExtendedNumber ExtendedNumber::ONE(){
    ExtendedNumber one(1, 0);
    return one;
}
