#include <gtest/gtest.h>
#include "number.hpp"
#include "polynomial.hpp"

class PolynomialTest : public ::testing::Test {
    protected:
    virtual void SetUp() {
        mpz_init_set_si(Number::CHARA, 31);
        mpz_init_set_si(Number::MCHARA, -31);
    }
};

TEST_F(PolynomialTest, InitDegreeEquals) {
    
    Number f0(3);
    Number f1(1);
    int d = 1;
    Polynomial<Number> f(d, f0, f1);
    EXPECT_EQ(d, f.deg);
}

TEST_F(PolynomialTest, SumEquals) {
    mpz_init_set_si(Number::CHARA, 31);

    // (4x^2 + x + 3) + (x + 2) = 4x^2 + 2x + 5.
    Number f0(3);
    Number f1(1);
    Number f2(4);
    int f_deg = 2;
    Polynomial<Number> f(f_deg);
    f.coeff[0] = f0; f.coeff[1] = f1; f.coeff[2] = f2;

    Number g0(2);
    Number g1(1);
    int g_deg = 1;
    Polynomial<Number> g(g_deg, g0, g1);

    Number h0(5);
    Number h1(2);
    Number h2(4);
    int h_deg = 2;
    Polynomial<Number> h(h_deg);
    h.coeff[0] = h0; h.coeff[1] = h1; h.coeff[2] = h2;

    EXPECT_EQ(h, f + g);
}

TEST_F(PolynomialTest, DifferenceEquals) {
    // (4x^2 + x + 3) - (x + 2) = 4x^2 + 1.
    Number f0(3);
    Number f1(1);
    Number f2(4);
    int f_deg = 2;
    Polynomial<Number> f(f_deg);
    f.coeff[0] = f0; f.coeff[1] = f1; f.coeff[2] = f2;

    Number g0(2);
    Number g1(1);
    int g_deg = 1;
    Polynomial<Number> g(g_deg, g0, g1);

    Number h0(1);
    Number h1(0);
    Number h2(4);
    int h_deg = 2;
    Polynomial<Number> h(h_deg);
    h.coeff[0] = h0; h.coeff[1] = h1; h.coeff[2] = h2;

    EXPECT_EQ(h, f - g);
}

TEST_F(PolynomialTest, ProductEquals) {
    // (4x^2 + x + 3) * (x + 2) = 4x^3 + 9x^2 + 5x + 6.
    Number f0(3);
    Number f1(1);
    Number f2(4);
    int f_deg = 2;
    Polynomial<Number> f(f_deg);
    f.coeff[0] = f0; f.coeff[1] = f1; f.coeff[2] = f2;

    Number g0(2);
    Number g1(1);
    int g_deg = 1;
    Polynomial<Number> g(g_deg, g0, g1);

    Number h0(6);
    Number h1(5);
    Number h2(9);
    Number h3(4);
    int h_deg = 3;
    Polynomial<Number> h(h_deg);
    h.coeff[0] = h0; h.coeff[1] = h1; h.coeff[2] = h2; h.coeff[3] = h3;

    EXPECT_EQ(h, f * g);
}

TEST_F(PolynomialTest, QuotientEquals) {
    // 4x^2 + x + 3 = (x + 2) * (4x - 7) + 17.
    Number f0(3);
    Number f1(1);
    Number f2(4);
    int f_deg = 2;
    Polynomial<Number> f(f_deg);
    f.coeff[0] = f0; f.coeff[1] = f1; f.coeff[2] = f2;

    Number g0(2);
    Number g1(1);
    int g_deg = 1;
    Polynomial<Number> g(g_deg, g0, g1);

    Number q0(-7);
    Number q1(4);
    int q_deg = 1;
    Polynomial<Number> q(q_deg, q0, q1);

    EXPECT_EQ(1, (f / g).deg);
    EXPECT_EQ(q, f / g);
}

TEST_F(PolynomialTest, RemainderEquals) {
    // 4x^2 + x + 3 = (x + 2) * (4x - 7) + 17.
    Number f0(3);
    Number f1(1);
    Number f2(4);
    int f_deg = 2;
    Polynomial<Number> f(f_deg);
    f.coeff[0] = f0; f.coeff[1] = f1; f.coeff[2] = f2;

    Number g0(2);
    Number g1(1);
    int g_deg = 1;
    Polynomial<Number> g(g_deg, g0, g1);

    Number r0(17);
    int r_deg = 0;
    Polynomial<Number> r(r_deg, r0);

    EXPECT_EQ(0, (f % g).deg);
    EXPECT_EQ(r, f % g);
}

TEST_F(PolynomialTest, DerivativeEquals) {
    // (4x^2 + x + 3)' = 8x + 1.
    Number f0(3);
    Number f1(1);
    Number f2(4);
    int f_deg = 2;
    Polynomial<Number> f(f_deg);
    f.coeff[0] = f0; f.coeff[1] = f1; f.coeff[2] = f2;

    Number g0(1);
    Number g1(8);
    int g_deg = 1;
    Polynomial<Number> g(g_deg, g0, g1);

    EXPECT_EQ(g, f.derivative());
}
