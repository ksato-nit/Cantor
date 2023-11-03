#include <gtest/gtest.h>
#include "number.hpp"
#include "polynomial.hpp"
#include "mumford.hpp"

TEST(PolynomialTest, SumOfDeg5MonicEquals) {
    int fc[6] = {-1, 3, 6, -2, -3, 1};
    int hc[3] = {0, 1, 1};
    int u1c[3] = {5, 25, 1};
    int v1c[2] = {-30, -27};
    int u2c[3] = {2, 11, 1};
    int v2c[2] = {-14, -21};
    int u12c[3] = {11, 4, 1};
    int v12c[2] = {-25, -8};

    Polynomial f(5, fc);
    Polynomial h(2, hc);
    Polynomial u1(2, u1c);
    Polynomial v1(1, v1c);
    Polynomial u2(2, u2c);
    Polynomial v2(1, v2c);
    Polynomial u12(2, u12c);
    Polynomial v12(1, v12c);

    Mumford D1(f, h, u1, v1);
    Mumford D2(f, h, u2, v2);
    Mumford D12(f, h, u12, v12);

    Mumford Sum = D1 + D2;

    EXPECT_EQ(D12.u, Sum.u);
    EXPECT_EQ(D12.v, Sum.v);
}

TEST(PolynomialTest, SumOfDeg5MonicDegeneratedEquals) {
    int fc[6] = {-1, 3, 6, -2, -3, 1};
    int hc[3] = {0, 1, 1};
    int u1c[3] = {5, 25, 1};
    int v1c[2] = {-30, -27};
    int u2c[3] = {17, 15, 1};
    int v2c[2] = {-25, -26};
    int u12c[2] = {13, 1};
    int v12c[1] = {0};

    Polynomial f(5, fc);
    Polynomial h(2, hc);
    Polynomial u1(2, u1c);
    Polynomial v1(1, v1c);
    Polynomial u2(2, u2c);
    Polynomial v2(1, v2c);
    Polynomial u12(1, u12c);
    Polynomial v12(0, v12c);

    Mumford D1(f, h, u1, v1);
    Mumford D2(f, h, u2, v2);
    Mumford D12(f, h, u12, v12);

    Mumford Sum = D1 + D2;

    EXPECT_EQ(D12.u, Sum.u);
    EXPECT_EQ(D12.v, Sum.v);
}

TEST(PolynomialTest, SumOfDeg6NonMonicEquals) {
    int fc[7] = {-1, 3, 6, -2, -3, 1, 1};
    int hc[1] = {0};
    int u1c[3] = {5, 25, 1};
    int v1c[2] = {-2, -23};
    int u2c[3] = {17, 17, 1};
    int v2c[2] = {-27, -16};
    int u12c[3] = {-14, 7, 1};
    int v12c[2] = {17, 21};    

    Polynomial f(6, fc);
    Polynomial h(0, hc);

    Polynomial u1(2, u1c);
    Polynomial v1(1, v1c);
    Polynomial u2(2, u2c);
    Polynomial v2(1, v2c);
    Polynomial u12(2, u12c);
    Polynomial v12(1, v12c);

    Mumford D1(f, h, u1, v1);
    Mumford D2(f, h, u2, v2);    

    Mumford D12(f, h, u12, v12);

    Mumford Sum = D1 + D2;

    EXPECT_EQ(D12.u, Sum.u);
    EXPECT_EQ(D12.v, Sum.v);
}
