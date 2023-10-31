#include "iostream"
#include "number.hpp"
#include "polynomial.hpp"
#include "mumford.hpp"
#include "mumford_projective.hpp"

int main(){
    int fc[6] = {-1, 3, 6, -2, -3, 1};
    int hc[1] = {0};
    int u1c[3] = {2, 28, 1};
    int v1c[2] = {-7, -22};
    int u2c[3] = {12, 24, 1};
    int v2c[2] = {-23, -8};

    Polynomial f(5, fc);
    Polynomial h(0, hc);

    Polynomial u1(2, u1c);
    Polynomial v1(1, v1c);
    Polynomial u2(2, u2c);
    Polynomial v2(1, v2c);

    Mumford D1(f, h, u1, v1);
    Mumford D2(f, h, u2, v2);

    std::cout << "D1:" << std::endl;
    D1.print();
    std::cout << "D2:" << std::endl;
    D2.print();

    std::cout << "D1 + D2:" << std::endl;
    Mumford sum1 = D1 + D2;
    sum1.print();

    return 0;
}
