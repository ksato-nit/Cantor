#include <gtest/gtest.h>
#include "number.hpp"
#include "polynomial.hpp"

TEST(PolynomialTest, InitDegreeEquals) {
    Number::CHARA = 31;
    Number::MCHARA = -31;
    
    Number f0(3);
    Number f1(1);
    int d = 1;
    Polynomial f(d, f0, f1);
    EXPECT_EQ(d, f.deg);
}

TEST(PolynomialTest, SumEquals) {
    Number::CHARA = 31;
    Number::MCHARA = -31;

    // (4x^2 + x + 3) + (x + 2) = 4x^2 + 2x + 5.
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

    Number h0(5);
    Number h1(2);
    Number h2(4);
    int h_deg = 2;
    Polynomial h(h_deg);
    h.coeff[0] = h0; h.coeff[1] = h1; h.coeff[2] = h2;

    EXPECT_EQ(h, f + g);
}

TEST(PolynomialTest, DifferenceEquals) {
    Number::CHARA = 31;
    Number::MCHARA = -31;

    // (4x^2 + x + 3) - (x + 2) = 4x^2 + 1.
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

    Number h0(1);
    Number h1(0);
    Number h2(4);
    int h_deg = 2;
    Polynomial h(h_deg);
    h.coeff[0] = h0; h.coeff[1] = h1; h.coeff[2] = h2;

    EXPECT_EQ(h, f - g);
}

TEST(PolynomialTest, ProductEquals) {
    Number::CHARA = 31;
    Number::MCHARA = -31;

    // (4x^2 + x + 3) * (x + 2) = 4x^3 + 9x^2 + 5x + 6.
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

    Number h0(6);
    Number h1(5);
    Number h2(9);
    Number h3(4);
    int h_deg = 3;
    Polynomial h(h_deg);
    h.coeff[0] = h0; h.coeff[1] = h1; h.coeff[2] = h2; h.coeff[3] = h3;

    EXPECT_EQ(h, f * g);
}

TEST(PolynomialTest, QuotientEquals) {
    Number::CHARA = 31;
    Number::MCHARA = -31;

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

    EXPECT_EQ(1, (f / g).deg);
    EXPECT_EQ(q, f / g);
}

TEST(PolynomialTest, RemainderEquals) {
    Number::CHARA = 31;
    Number::MCHARA = -31;

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

    EXPECT_EQ(0, (f % g).deg);
    EXPECT_EQ(r, f % g);
}

TEST(PolynomialTest, DerivativeEquals) {
    Number::CHARA = 31;
    Number::MCHARA = -31;

    // (4x^2 + x + 3)' = 8x + 1.
    Number f0(3);
    Number f1(1);
    Number f2(4);
    int f_deg = 2;
    Polynomial f(f_deg);
    f.coeff[0] = f0; f.coeff[1] = f1; f.coeff[2] = f2;

    Number g0(1);
    Number g1(8);
    int g_deg = 1;
    Polynomial g(g_deg, g0, g1);

    EXPECT_EQ(g, f.derivative());
}
