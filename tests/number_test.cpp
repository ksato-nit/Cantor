#include <gtest/gtest.h>
#include "number.hpp"

TEST(NumberTest, ModValueEquals) {
    Number::CHARA = 31;
    Number::MCHARA = -31;
    Number x(1);
    Number y(1);
    Number C(31);
    y = y - C;
    EXPECT_EQ(x, y);
}

TEST(NumberTest, ModValueNotEquals) {
    Number::CHARA = 31;
    Number::MCHARA = -31;
    Number x(1);
    Number y(2);
    Number C(31);
    y = y - C;
    EXPECT_NE(x, y);
}

TEST(NumberTest, SumEquals) {
    Number::CHARA = 31;
    Number::MCHARA = -31;
    Number x(1);
    Number y(2);
    Number z(3);
    EXPECT_EQ(z, x + y);
}

TEST(NumberTest, LargeSumEquals) {
    Number::CHARA = 31;
    Number::MCHARA = -31;
    Number x(50);
    Number y(40);
    Number z(90);
    EXPECT_EQ(z, x + y);
}

TEST(NumberTest, DifferenceEquals) {
    Number::CHARA = 31;
    Number::MCHARA = -31;
    Number x(5);
    Number y(4);
    Number z(1);
    EXPECT_EQ(z, x - y);
}

TEST(NumberTest, LargeDifferenceEquals) {
    Number::CHARA = 31;
    Number::MCHARA = -31;
    Number x(200);
    Number y(40);
    Number z(160);
    EXPECT_EQ(z, x - y);
}

TEST(NumberTest, ProductEquals) {
    Number::CHARA = 31;
    Number::MCHARA = -31;
    Number x(2);
    Number y(3);
    Number z(6);
    EXPECT_EQ(z, x * y);
}

TEST(NumberTest, LargeProductEquals) {
    Number::CHARA = 31;
    Number::MCHARA = -31;
    Number x(50);
    Number y(40);
    Number z(2000);
    EXPECT_EQ(z, x * y);
}

TEST(NumberTest, ProductToIntEquals) {
    Number::CHARA = 31;
    Number::MCHARA = -31;
    Number x(9);
    int y = 2;
    Number z(18);
    EXPECT_EQ(z, x * y);
}

TEST(NumberTest, QuotientEquals) {
    Number::CHARA = 31;
    Number::MCHARA = -31;
    Number x(2);
    Number y(3);
    Number z = x / y;
    EXPECT_EQ(x, z * y);
}
