#include <gtest/gtest.h>
#include "number.hpp"
#include "polynomial.hpp"
#include "mumford_projective.hpp"

TEST(ProjectiveMumfordTest, SumEquals) {
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

    Polynomial u1_half = u1 * Number(2);
    Polynomial v1_half = v1 * Number(2);
    Polynomial u2_half = u2 * Number(3);
    Polynomial v2_half = v2 * Number(3);

    ProjectiveMumford D1P(f, h, u1_half.coeff[1], u1_half.coeff[0], v1_half.coeff[1], v1_half.coeff[0], Number(2));
    ProjectiveMumford D2P(f, h, u2_half.coeff[1], u2_half.coeff[0], v2_half.coeff[1], v2_half.coeff[0], Number(3));
    ProjectiveMumford Sum = D1P + D2P;

    ProjectiveMumford D12(f, h, u12.coeff[1], u12.coeff[0], v12.coeff[1], v12.coeff[0]);

    EXPECT_EQ(D12.U1, Sum.U1 / Sum.Z);
    EXPECT_EQ(D12.U0, Sum.U0 / Sum.Z);
    EXPECT_EQ(D12.V1, Sum.V1 / Sum.Z);
    EXPECT_EQ(D12.V0, Sum.V0 / Sum.Z);    
}