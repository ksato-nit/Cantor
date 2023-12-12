#include <gtest/gtest.h>
#include "number.hpp"
#include "polynomial.hpp"
#include "mumford.hpp"

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

    Mumford D1(f, h, u1, v1);
    Mumford D2(f, h, u2, v2);
    Mumford D12(f, h, u12, v12);

    Mumford Sum = D1 + D2;

    EXPECT_EQ(D12.u, Sum.u);
    EXPECT_EQ(D12.v, Sum.v);
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

    Mumford D1(f, h, u1, v1);
    Mumford D2(f, h, u2, v2);
    Mumford D12(f, h, u12, v12);

    Mumford Sum = D1 + D2;

    EXPECT_EQ(D12.u, Sum.u);
    EXPECT_EQ(D12.v, Sum.v);
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

    Mumford D1(f, h, u1, v1);
    Mumford D2(f, h, u2, v2);    

    Mumford D12(f, h, u12, v12);

    Mumford Sum = D1.CantorAdd(D2);

    EXPECT_EQ(D12.u, Sum.u);
    EXPECT_EQ(D12.v, Sum.v);
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

    Mumford D1(f, h, u1, v1);
    Mumford D2(f, h, u2, v2);    

    Mumford D12(f, h, u12, v12);

    Mumford Sum = D1.LangeAdd(D2);

    EXPECT_EQ(D12.u, Sum.u);
    EXPECT_EQ(D12.v, Sum.v);
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

    Mumford D1(f, h, u1, v1);
    Mumford D2(f, h, u2, v2);    

    Mumford D12(f, h, u12, v12);

    Mumford Sum = D1.CostelloAdd(D2);

    EXPECT_EQ(D12.u, Sum.u);
    EXPECT_EQ(D12.v, Sum.v);
}

TEST(MumfordTest, DoublingDeg6NonMonicEquals) {
    int fc[7] = {-1, 3, 6, -2, -3, 1, 1};
    int hc[1] = {0};
    int uc[3] = {6, 24, 1};
    int vc[2] = {-24, -1};
    int udc[3] = {-28, -14, 1};
    int vdc[2] = {-10, -13};    

    Polynomial f(6, fc);
    Polynomial h(0, hc);

    Polynomial u(2, uc);
    Polynomial v(1, vc);
    Polynomial ud(2, udc);
    Polynomial vd(1, vdc);

    Mumford D1(f, h, u, v);
    Mumford Ret(f, h, ud, vd);

    Mumford Sum = D1.doubling();

    EXPECT_EQ(Ret.u, Sum.u);
    EXPECT_EQ(Ret.v, Sum.v);
}
