#include <gtest/gtest.h>
#include "number.hpp"

TEST(NumberTest, SumEquals) {
    Number x(1);
    Number y(2);
    Number z(3);
    EXPECT_EQ(z, x + y);
}

TEST(NumberTest, DifferenceEquals) {
    Number x(5);
    Number y(4);
    Number z(1);
    EXPECT_EQ(z, x - y);
}
