#include "extended_number.hpp"

mpz_t Number::CHARA;

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

bool ExtendedNumber::isZero(){
    mpz_t zero;
    mpz_init_set_ui(zero, 0);
    int comp1 = mpz_comp(this->re, zero);
    int comp2 = mpz_comp(this->re, zero);
    return (comp1 == 0) && (comp2 == 0);
}

ExtendedNumber ExtendedNumber::ZERO(){
    ExtendedNumber one(0, 0);
    return zero;
}

ExtendedNumber ExtendedNumber::ONE(){
    ExtendedNumber one(1, 0);
    return one;
}
