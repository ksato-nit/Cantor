#include <gtest/gtest.h>
#include "polynomial.hpp"
#include "extended_number.hpp"
#include <gmpxx.h>

TEST(ExtendedNumberTest, ModValueEquals) {
    mpz_init_set_si(Number::CHARA, 31);
    mpz_init_set_si(Number::MCHARA, -31);

    Polynomial f(2);
    f.coeff[0] = Number::ONE(); f.coeff[1] = Number::ONE(); f.coeff[2] = Number::ONE();
    ExtendedNumber x(f, Number::ONE(), Number::ZERO());
    ExtendedNumber y(f, Number::ONE(), Number::ZERO()); // todo: なんとかする
    EXPECT_EQ(x, y);
}

TEST(ExtendedNumberTest, SumEquals) {
    mpz_init_set_si(Number::CHARA, 31);
    mpz_init_set_si(Number::MCHARA, -31);

    Polynomial f(2);
    f.coeff[0] = Number::ONE(); f.coeff[1] = Number::ONE(); f.coeff[2] = Number::ONE();

    Polynomial g1(1, Number(1), Number(1));
    Polynomial g2(1, Number(2), Number(0));
    Polynomial g3(1, Number(3), Number(1));
    ExtendedNumber x(f, g1);
    ExtendedNumber y(f, g2);
    ExtendedNumber z(f, g3);
    EXPECT_EQ(z, x + y);
}

TEST(ExtendedNumberTest, ProductEquals) {
    mpz_init_set_si(Number::CHARA, 31);
    mpz_init_set_si(Number::MCHARA, -31);

    Polynomial f(2);
    f.coeff[0] = Number::ONE(); f.coeff[1] = Number::ZERO(); f.coeff[2] = Number::ONE();

    Polynomial g1(1, Number(1), Number(1));
    Polynomial g2(1, Number(0), Number(2));
    Polynomial g3(1, Number(-2), Number(2));
    ExtendedNumber x(f, g1);
    ExtendedNumber y(f, g2);
    ExtendedNumber z(f, g3);
    EXPECT_EQ(z, x * y);
}

TEST(ExtendedNumberTest, InverseEquals) {
    mpz_init_set_si(Number::CHARA, 31);
    mpz_init_set_si(Number::MCHARA, -31);

    Polynomial f(2);
    f.coeff[0] = Number::ONE(); f.coeff[1] = Number::ZERO(); f.coeff[2] = Number::ONE();

    Polynomial g1(1, Number::ONE(), Number::ONE());
    Polynomial g2(0, Number::ONE());
    ExtendedNumber x(f, g1);
    ExtendedNumber y(f, g2);
    EXPECT_EQ(y, x.inv() * x);
}
