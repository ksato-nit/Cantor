#include <gtest/gtest.h>
#include "number.hpp"
#include "polynomial.hpp"
#include "iostream"

TEST(PolynomialTest, InitDegreeEquals) {
    Number f0(3);
    Number f1(1);
    int d = 1;
    Polynomial f(d, f0, f1);
    EXPECT_EQ(d, f.deg);
}

TEST(PolynomialTest, QuotientEquals) {
    // 4x^2 + x + 3 = (x + 2) * (4x - 7) + 17.
    Number f0(3);
    Number f1(1);
    Number f2(4);
    int f_deg = 2;
    Polynomial f(f_deg);
    f.coeff[0] = f0; f.coeff[1] = f1; f.coeff[2] = f2;

    Number g0(2);
    Number g1(1);
    int g_deg = 1;
    Polynomial g(g_deg, g0, g1);

    Number q0(-7);
    Number q1(4);
    int q_deg = 1;
    Polynomial q(q_deg, q0, q1);

    EXPECT_EQ(q, f / g);
}

TEST(PolynomialTest, RemainderEquals) {
    // 4x^2 + x + 3 = (x + 2) * (4x - 7) + 17.
    Number f0(3);
    Number f1(1);
    Number f2(4);
    int f_deg = 2;
    Polynomial f(f_deg);
    f.coeff[0] = f0; f.coeff[1] = f1; f.coeff[2] = f2;

    Number g0(2);
    Number g1(1);
    int g_deg = 1;
    Polynomial g(g_deg, g0, g1);

    Number r0(17);
    int r_deg = 0;
    Polynomial r(r_deg, r0);

    EXPECT_EQ(r, f % g);
}
