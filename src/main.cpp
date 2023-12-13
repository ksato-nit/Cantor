#include <iostream>
#include <chrono>
#include "number.hpp"
#include "polynomial.hpp"
#include "mumford.hpp"
#include "mumford_projective.hpp"
#include "mumford_weighted_projective.hpp"

int main(){
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

    Mumford D1(f, h, u1, v1);
    Mumford sum1 = D1.doubling();
    sum1.print();

    Polynomial u1_half = u1 * Number(2);
    Polynomial v1_half = v1 * Number(2);

    ProjectiveMumford D1P(f, h, u1_half.coeff[1], u1_half.coeff[0], v1_half.coeff[1], v1_half.coeff[0], Number(2));
    ProjectiveMumford sum2 = D1P.doubling();
    sum2.print();

    return 0;
}
