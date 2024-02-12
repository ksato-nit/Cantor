#include <gtest/gtest.h>
#include "number.hpp"

class NumberTest : public ::testing::Test {
    protected:
    virtual void SetUp() {
        mpz_init_set_si(Number::CHARA, 31);
        mpz_init_set_si(Number::MCHARA, -31);
    }
};

TEST_F(NumberTest, ModValueEquals) {
    Number x(1);
    Number y(1);
    Number C(31);
    y = y - C;
    EXPECT_EQ(x, y);
}

TEST_F(NumberTest, ModValueNotEquals) {
    Number x(1);
    Number y(2);
    Number C(31);
    y = y - C;
    EXPECT_NE(x, y);
}

TEST_F(NumberTest, SumEquals) {
    Number x(1);
    Number y(2);
    Number z(3);
    EXPECT_EQ(z, x + y);
}

TEST_F(NumberTest, LargeSumEquals) {
    Number x(50);
    Number y(40);
    Number z(90);
    EXPECT_EQ(z, x + y);
}

TEST_F(NumberTest, DifferenceEquals) {
    Number x(5);
    Number y(4);
    Number z(1);
    EXPECT_EQ(z, x - y);
}

TEST_F(NumberTest, LargeDifferenceEquals) {
    Number x(200);
    Number y(40);
    Number z(160);
    EXPECT_EQ(z, x - y);
}

TEST_F(NumberTest, ProductEquals) {
    Number x(2);
    Number y(3);
    Number z(6);
    EXPECT_EQ(z, x * y);
}

TEST_F(NumberTest, LargeProductEquals) {
    Number x(50);
    Number y(40);
    Number z(2000);
    EXPECT_EQ(z, x * y);
}

TEST_F(NumberTest, ProductToIntEquals) {
    Number x(9);
    int y = 2;
    Number z(18);
    EXPECT_EQ(z, x * y);
}

TEST_F(NumberTest, QuotientEquals) {
    Number x(2);
    Number y(3);
    Number z = x / y;
    EXPECT_EQ(x, z * y);
}
