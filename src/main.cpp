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

    Mumford sum = D1.CostelloDoubling();
    sum.print();

    return 0;
}
