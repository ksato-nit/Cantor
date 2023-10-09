#include "iostream"
#include "number.hpp"
#include "polynomial.hpp"
#include "mumford.hpp"

int main(){
    int fc[6] = {-1, 3, 6, -2, -3, 1};
    int hc[3] = {0, 1, 1};
    int u1c[2] = {30, 1};
    int v1c[1] = {-26};
    int u2c[3] = {30, 20, 1};
    int v2c[2] = {-29, -21};

    Polynomial f(5, fc);
    Polynomial h(2, hc);
    Polynomial u1(1, u1c);
    Polynomial v1(0, v1c);
    Polynomial u2(2, u2c);
    Polynomial v2(1, v2c);

    Mumford D1(f, h, u1, v1);
    Mumford D2(f, h, u2, v2);

    Mumford sum = D1 + D2;
    sum.print();

    return 0;
}
