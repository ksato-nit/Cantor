#include <gtest/gtest.h>
#include "number.hpp"
#include "polynomial.hpp"
#include "mumford.hpp"

// fixme: Mumford の設計変更に伴い，回りくどくなっている．

TEST(MumfordTest, SumOfDeg5MonicEquals) {
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

    Mumford D1(f, h, u1.coeff[1], u1.coeff[0], v1.coeff[1], v1.coeff[0]);
    Mumford D2(f, h, u2.coeff[1], u2.coeff[0], v2.coeff[1], v2.coeff[0]);
    Mumford D12(f, h, u12.coeff[1], u12.coeff[0], v12.coeff[1], v12.coeff[0]);

    Mumford Sum = D1 + D2;

    EXPECT_EQ(D12.u1, Sum.u1);
    EXPECT_EQ(D12.u0, Sum.u0);
    EXPECT_EQ(D12.v1, Sum.v1);
    EXPECT_EQ(D12.v0, Sum.v0);
}

TEST(MumfordTest, SumOfDeg5MonicDegeneratedEquals) {
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

    Mumford D1(f, h, u1.coeff[1], u1.coeff[0], v1.coeff[1], v1.coeff[0]);
    Mumford D2(f, h, u2.coeff[1], u2.coeff[0], v2.coeff[1], v2.coeff[0]);
    Mumford D12(f, h, u12.coeff[1], u12.coeff[0], v12.coeff[1], v12.coeff[0]);

    Mumford Sum = D1 + D2;

    EXPECT_EQ(D12.u1, Sum.u1);
    EXPECT_EQ(D12.u0, Sum.u0);
    EXPECT_EQ(D12.v1, Sum.v1);
    EXPECT_EQ(D12.v0, Sum.v0);
}

TEST(MumfordTest, CantorSumOfDeg6NonMonicEquals) {
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

    Mumford D1(f, h, u1.coeff[1], u1.coeff[0], v1.coeff[1], v1.coeff[0]);
    Mumford D2(f, h, u2.coeff[1], u2.coeff[0], v2.coeff[1], v2.coeff[0]);
    Mumford D12(f, h, u12.coeff[1], u12.coeff[0], v12.coeff[1], v12.coeff[0]);

    Mumford Sum = D1.CantorAdd(D2);

    EXPECT_EQ(D12.u1, Sum.u1);
    EXPECT_EQ(D12.u0, Sum.u0);
    EXPECT_EQ(D12.v1, Sum.v1);
    EXPECT_EQ(D12.v0, Sum.v0);
}

TEST(MumfordTest, LangeSumOfDeg6NonMonicEquals) {
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

    Mumford D1(f, h, u1.coeff[1], u1.coeff[0], v1.coeff[1], v1.coeff[0]);
    Mumford D2(f, h, u2.coeff[1], u2.coeff[0], v2.coeff[1], v2.coeff[0]);
    Mumford D12(f, h, u12.coeff[1], u12.coeff[0], v12.coeff[1], v12.coeff[0]);

    Mumford Sum = D1.LangeAdd(D2);

    EXPECT_EQ(D12.u1, Sum.u1);
    EXPECT_EQ(D12.u0, Sum.u0);
    EXPECT_EQ(D12.v1, Sum.v1);
    EXPECT_EQ(D12.v0, Sum.v0);
}

TEST(MumfordTest, CostelloSumOfDeg6NonMonicEquals) {
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

    Mumford D1(f, h, u1.coeff[1], u1.coeff[0], v1.coeff[1], v1.coeff[0]);
    Mumford D2(f, h, u2.coeff[1], u2.coeff[0], v2.coeff[1], v2.coeff[0]);
    Mumford D12(f, h, u12.coeff[1], u12.coeff[0], v12.coeff[1], v12.coeff[0]);

    Mumford Sum = D1.CostelloAdd(D2);

    EXPECT_EQ(D12.u1, Sum.u1);
    EXPECT_EQ(D12.u0, Sum.u0);
    EXPECT_EQ(D12.v1, Sum.v1);
    EXPECT_EQ(D12.v0, Sum.v0);
}

TEST(MumfordTest, DoublingDeg6NonMonicEquals) {
    int fc[7] = {-1, 3, 6, -2, -3, 1, 1};
    int hc[1] = {0};
    int uc[3] = {6, 24, 1};
    int vc[2] = {-24, -1};
    int udc[3] = {-28, -14, 1};
    int vdc[2] = {10, 13};    

    Polynomial f(6, fc);
    Polynomial h(0, hc);

    Polynomial u(2, uc);
    Polynomial v(1, vc);
    Polynomial ud(2, udc);
    Polynomial vd(1, vdc);

    Mumford D1(f, h, u.coeff[1], u.coeff[0], v.coeff[1], v.coeff[0]);
    Mumford D12(f, h, ud.coeff[1], ud.coeff[0], vd.coeff[1], vd.coeff[0]);

    Mumford Sum = D1.LangeDoubling();

    EXPECT_EQ(D12.u1, Sum.u1);
    EXPECT_EQ(D12.u0, Sum.u0);
    EXPECT_EQ(D12.v1, Sum.v1);
    EXPECT_EQ(D12.v0, Sum.v0);
}

TEST(MumfordTest, ScalarMultiplicationEquals) {
    int fc[7] = {-1, 3, 6, -2, -3, 1, 1};
    int hc[1] = {0};

    Polynomial f(6, fc);
    Polynomial h(0, hc);

    Polynomial u1(2);
    Polynomial v1(1);

    u1.coeff[2].value.set_str("1", 10);
    u1.coeff[1].value.set_str("25", 10);
    u1.coeff[0].value.set_str("5", 10);

    v1.coeff[1].value.set_str("-23", 10);
    v1.coeff[0].value.set_str("-2", 10);

    Mumford D1(f, h, u1.coeff[1], u1.coeff[0], v1.coeff[1], v1.coeff[0]);

    mpz_class k = 200;
    Mumford Dk = D1 * k;

    Polynomial u2(2);
    Polynomial v2(1);

    u2.coeff[2].value.set_str("1", 10);
    u2.coeff[1].value.set_str("11", 10);
    u2.coeff[0].value.set_str("28", 10);

    v2.coeff[1].value.set_str("-29", 10);
    v2.coeff[0].value.set_str("-21", 10);

    EXPECT_EQ(Dk.u1, u2.coeff[1]);
    EXPECT_EQ(Dk.u0, u2.coeff[0]);
    EXPECT_EQ(Dk.v1, v2.coeff[1]);
    EXPECT_EQ(Dk.v0, v2.coeff[0]);
}
