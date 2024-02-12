#include <gtest/gtest.h>
#include "number.hpp"
#include "polynomial.hpp"
#include "mumford_projective.hpp"

TEST(ProjectiveMumfordTest, LangeSumEquals) {
    mpz_init_set_si(Number::CHARA, 31);
    mpz_init_set_si(Number::MCHARA, -31);

    int fc[7] = {-1, 3, 6, -2, -3, 1, 1};
    int hc[1] = {0};
    int u1c[3] = {5, 25, 1};
    int v1c[2] = {-2, -23};
    int u2c[3] = {17, 17, 1};
    int v2c[2] = {-27, -16};
    int u12c[3] = {-14, 7, 1};
    int v12c[2] = {17, 21};

    Polynomial<Number> f(6, fc);
    Polynomial<Number> h(0, hc);
    Polynomial<Number> u1(2, u1c);
    Polynomial<Number> v1(1, v1c);
    Polynomial<Number> u2(2, u2c);
    Polynomial<Number> v2(1, v2c);
    Polynomial<Number> u12(2, u12c);
    Polynomial<Number> v12(1, v12c);

    Polynomial<Number> u1_half = u1 * Number(2);
    Polynomial<Number> v1_half = v1 * Number(2);
    Polynomial<Number> u2_half = u2 * Number(3);
    Polynomial<Number> v2_half = v2 * Number(3);

    ProjectiveMumford<Number> D1P(f, h, u1_half.coeff[1], u1_half.coeff[0], v1_half.coeff[1], v1_half.coeff[0], Number(2));
    ProjectiveMumford<Number> D2P(f, h, u2_half.coeff[1], u2_half.coeff[0], v2_half.coeff[1], v2_half.coeff[0], Number(3));
    ProjectiveMumford<Number> Sum = D1P.LangeAdd(D2P);

    ProjectiveMumford<Number> D12(f, h, u12.coeff[1], u12.coeff[0], v12.coeff[1], v12.coeff[0]);

    EXPECT_EQ(D12.U1, Sum.U1 / Sum.Z);
    EXPECT_EQ(D12.U0, Sum.U0 / Sum.Z);
    EXPECT_EQ(D12.V1, Sum.V1 / Sum.Z);
    EXPECT_EQ(D12.V0, Sum.V0 / Sum.Z);    
}

TEST(ProjectiveMumfordTest, CostelloSumEquals) {
    mpz_init_set_si(Number::CHARA, 31);
    mpz_init_set_si(Number::MCHARA, -31);

    int fc[7] = {-1, 3, 6, -2, -3, 1, 1};
    int hc[1] = {0};
    int u1c[3] = {5, 25, 1};
    int v1c[2] = {-2, -23};
    int u2c[3] = {17, 17, 1};
    int v2c[2] = {-27, -16};
    int u12c[3] = {-14, 7, 1};
    int v12c[2] = {17, 21};

    Polynomial<Number> f(6, fc);
    Polynomial<Number> h(0, hc);
    Polynomial<Number> u1(2, u1c);
    Polynomial<Number> v1(1, v1c);
    Polynomial<Number> u2(2, u2c);
    Polynomial<Number> v2(1, v2c);
    Polynomial<Number> u12(2, u12c);
    Polynomial<Number> v12(1, v12c);

    Polynomial<Number> u1_half = u1 * Number(2);
    Polynomial<Number> v1_half = v1 * Number(2);
    Polynomial<Number> u2_half = u2 * Number(3);
    Polynomial<Number> v2_half = v2 * Number(3);

    ProjectiveMumford<Number> D1P(f, h, u1_half.coeff[1], u1_half.coeff[0], v1_half.coeff[1], v1_half.coeff[0], Number(2));
    ProjectiveMumford<Number> D2P(f, h, u2_half.coeff[1], u2_half.coeff[0], v2_half.coeff[1], v2_half.coeff[0], Number(3));
    ProjectiveMumford<Number> Sum = D1P.CostelloAdd(D2P);

    ProjectiveMumford<Number> D12(f, h, u12.coeff[1], u12.coeff[0], v12.coeff[1], v12.coeff[0]);

    EXPECT_EQ(D12.U1, Sum.U1 / Sum.Z);
    EXPECT_EQ(D12.U0, Sum.U0 / Sum.Z);
    EXPECT_EQ(D12.V1, Sum.V1 / Sum.Z);
    EXPECT_EQ(D12.V0, Sum.V0 / Sum.Z);    
}

TEST(ProjectiveMumfordTest, DoublingEquals) {
    mpz_init_set_si(Number::CHARA, 31);
    mpz_init_set_si(Number::MCHARA, -31);
    
    int fc[7] = {-1, 3, 6, -2, -3, 1, 1};
    int hc[1] = {0};
    int uc[3] = {6, 24, 1};
    int vc[2] = {-24, -1};
    int udc[3] = {-28, -14, 1};
    int vdc[2] = {10, 13};    

    Polynomial<Number> f(6, fc);
    Polynomial<Number> h(0, hc);

    Polynomial<Number> u(2, uc);
    Polynomial<Number> v(1, vc);
    Polynomial<Number> ud(2, udc);
    Polynomial<Number> vd(1, vdc);

    Polynomial<Number> u_half = u * Number(2);
    Polynomial<Number> v_half = v * Number(2);

    ProjectiveMumford<Number> DP(f, h, u_half.coeff[1], u_half.coeff[0], v_half.coeff[1], v_half.coeff[0], Number(2));
    ProjectiveMumford<Number> Sum = DP.LangeDoubling();

    ProjectiveMumford<Number> D12(f, h, ud.coeff[1], ud.coeff[0], vd.coeff[1], vd.coeff[0]);

    EXPECT_EQ(D12.U1, Sum.U1 / Sum.Z);
    EXPECT_EQ(D12.U0, Sum.U0 / Sum.Z);
    EXPECT_EQ(D12.V1, Sum.V1 / Sum.Z);
    EXPECT_EQ(D12.V0, Sum.V0 / Sum.Z); 
}
