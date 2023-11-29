#include <gtest/gtest.h>
#include "polynomial.hpp"
#include "extended_number.hpp"
#include <gmpxx.h>

TEST(ExtehdedNumberTest, ModValueEquals) {
    Polynomial f(2);
    f.coeff[0] = Number::ONE(); f.coeff[1] = Number::ONE(); f.coeff[2] = Number::ONE();
    ExtendedNumber x(f, Number::ONE(), Number::ZERO());
    ExtendedNumber y(f, Number::ONE() - x.CHARA, Number::ZERO());
    EXPECT_EQ(x, y);
}

TEST(ExtehdedNumberTest, SumEquals) {
    Polynomial f(2);
    f.coeff[0] = Number::ONE(); f.coeff[1] = Number::ONE(); f.coeff[2] = Number::ONE();

    Polynomial g1(1, Number(mpz_class(1)), Number(mpz_class(1)));
    Polynomial g2(1, Number(mpz_class(2)), Number(mpz_class(0)));
    Polynomial g3(1, Number(mpz_class(3)), Number(mpz_class(1)));
    ExtendedNumber x(f, g1);
    ExtendedNumber y(f, g2);
    ExtendedNumber z(f, g3);
    EXPECT_EQ(z, x + y);
}

TEST(ExtehdedNumberTest, ProductEquals) {
    Polynomial f(2);
    f.coeff[0] = Number::ONE(); f.coeff[1] = Number::ZERO(); f.coeff[2] = Number::ONE();

    Polynomial g1(1, Number(mpz_class(1)), Number(mpz_class(1)));
    Polynomial g2(1, Number(mpz_class(0)), Number(mpz_class(2)));
    Polynomial g3(1, Number(mpz_class(-2)), Number(mpz_class(2)));
    ExtendedNumber x(f, g1);
    ExtendedNumber y(f, g2);
    ExtendedNumber z(f, g3);
    EXPECT_EQ(z, x * y);
}

TEST(ExtehdedNumberTest, InverseEquals) {
    Polynomial f(2);
    f.coeff[0] = Number::ONE(); f.coeff[1] = Number::ZERO(); f.coeff[2] = Number::ONE();

    Polynomial g1(1, Number(mpz_class(1)), Number(mpz_class(1)));
    Polynomial g2(0, Number::ONE());
    ExtendedNumber x(f, g1);
    ExtendedNumber y(f, g2);
    EXPECT_EQ(y, x.inv() * x);
}
