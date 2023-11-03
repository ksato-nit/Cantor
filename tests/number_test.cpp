#include <gtest/gtest.h>
#include "number.hpp"

TEST(NumberTest, ModValueEquals) {
    Number x(1);
    Number y(1 - Number::CHARA);
    EXPECT_EQ(x, y);
}

TEST(NumberTest, SumEquals) {
    Number x(1);
    Number y(2);
    Number z(3);
    EXPECT_EQ(z, x + y);
}

TEST(NumberTest, LargeSumEquals) {
    Number x(50);
    Number y(40);
    Number z(90);
    EXPECT_EQ(z, x + y);
}

TEST(NumberTest, DifferenceEquals) {
    Number x(5);
    Number y(4);
    Number z(1);
    EXPECT_EQ(z, x - y);
}

TEST(NumberTest, LargeDifferenceEquals) {
    Number x(200);
    Number y(40);
    Number z(160);
    EXPECT_EQ(z, x - y);
}

TEST(NumberTest, ProductEquals) {
    Number x(2);
    Number y(3);
    Number z(6);
    EXPECT_EQ(z, x * y);
}

TEST(NumberTest, LargeProductEquals) {
    Number x(50);
    Number y(40);
    Number z(2000);
    EXPECT_EQ(z, x * y);
}

TEST(NumberTest, ProductToIntEquals) {
    Number x(9);
    int y = 2;
    Number z(18);
    EXPECT_EQ(z, x * y);
}

TEST(NumberTest, QuotientEquals) {
    Number x(2);
    Number y(3);
    Number z = x / y;
    EXPECT_EQ(x, z * y);
}
