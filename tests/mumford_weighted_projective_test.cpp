#include <gtest/gtest.h>
#include "number.hpp"
#include "polynomial.hpp"
#include "mumford_weighted_projective.hpp"

TEST(WeightedProjectiveMumfordTest, SumEquals) {
    Number::CHARA.set_str("31", 10);
    
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

    Number Z11 = Number(2);
    Number Z12 = Number(3);
    Number Z21 = Number(3);
    Number Z22 = Number(5);
    Polynomial u1_half = u1 * Z11 * Z11;
    Polynomial v1_half = v1 * Z11 * Z11 * Z11 * Z12;
    Polynomial u2_half = u2 * Z21 * Z21;
    Polynomial v2_half = v2 * Z21 * Z21 * Z21 * Z22;

    WeightedProjectiveMumford D1WP(f, h, u1_half.coeff[1], u1_half.coeff[0], v1_half.coeff[1], v1_half.coeff[0], Z11, Z12);
    WeightedProjectiveMumford D2WP(f, h, u2_half.coeff[1], u2_half.coeff[0], v2_half.coeff[1], v2_half.coeff[0], Z21, Z22);
    WeightedProjectiveMumford Sum = D1WP + D2WP;

    WeightedProjectiveMumford D12(f, h, u12.coeff[1], u12.coeff[0], v12.coeff[1], v12.coeff[0]);

    EXPECT_EQ(D12.U1, Sum.U1 / (Sum.Z1 * Sum.Z1));
    EXPECT_EQ(D12.U0, Sum.U0 / (Sum.Z1 * Sum.Z1));
    EXPECT_EQ(D12.V1, Sum.V1 / (Sum.Z1 * Sum.Z1 * Sum.Z1 * Sum.Z2));
    EXPECT_EQ(D12.V0, Sum.V0 / (Sum.Z1 * Sum.Z1 * Sum.Z1 * Sum.Z2));    
}
