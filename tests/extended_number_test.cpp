#include <gtest/gtest.h>
#include "polynomial.hpp"
#include "extended_number.hpp"
#include <gmpxx.h>

TEST(ExtendedNumberTest, SumEquals) {
    mpz_init_set_si(Number::CHARA, 31);
    mpz_init_set_si(Number::MCHARA, -31);

    ExtendedNumber x(Number(1), Number(2));
    ExtendedNumber y(Number(2), Number(3));
    ExtendedNumber z(Number(3), Number(5));
    EXPECT_EQ(z, x + y);
}


TEST(ExtendedNumberTest, ProductEquals) {
    mpz_init_set_si(Number::CHARA, 31);
    mpz_init_set_si(Number::MCHARA, -31);

    ExtendedNumber x(Number(1), Number(2));
    ExtendedNumber y(Number(2), Number(3));
    ExtendedNumber z(Number(-4), Number(7));
    EXPECT_EQ(z, x * y);
}

TEST(ExtendedNumberTest, InverseEquals) {
    mpz_init_set_si(Number::CHARA, 31);
    mpz_init_set_si(Number::MCHARA, -31);

    ExtendedNumber x(Number(1), Number(2));
    EXPECT_EQ(ExtendedNumber::ONE(), x.inv() * x);
}
