#include "iostream"
#include "number.hpp"
#include "polynomial.hpp"
#include "mumford.hpp"
#include "mumford_projective.hpp"

int main(){
    int fc[6] = {-1, 3, 6, -2, -3, 1};
    int hc[3] = {0, 1, 1};
    int u1c[3] = {4, 26, 1};
    int v1c[2] = {-2, -24};
    int u2c[3] = {30, 20, 1};
    int v2c[2] = {-29, -21};
    int u3c[2] = {30, 1};
    int v3c[1] = {-26};

    Polynomial f(5, fc);
    Polynomial h(2, hc);

    Polynomial u1(2, u1c);
    Polynomial v1(1, v1c);
    Polynomial u2(2, u2c);
    Polynomial v2(1, v2c);
    Polynomial u3(1, u3c);
    Polynomial v3(0, v3c);

    Mumford D1(f, h, u1, v1);
    Mumford D2(f, h, u2, v2);
    Mumford D3(f, h, u3, v3);

    std::cout << "D1:" << std::endl;
    D1.print();
    std::cout << "-D1:" << std::endl;
    D1.inv().print();
    std::cout << "D2:" << std::endl;
    D2.print();
    std::cout << "D3:" << std::endl;
    D3.print();        

    std::cout << "D1 + D2:" << std::endl;
    Mumford sum1 = D1 + D2;
    sum1.print();

    ProjectiveMumford D1P(f, h, u1.coeff[1], u1.coeff[0], v1.coeff[1], v1.coeff[0]);
    ProjectiveMumford D2P(f, h, u2.coeff[1], u2.coeff[0], v2.coeff[1], v2.coeff[0]);
    ProjectiveMumford D3P(f, h, u3.coeff[1], u3.coeff[0], v3.coeff[1], v3.coeff[0]);
    std::cout << "D1P:" << std::endl;
    D1P.print();
    std::cout << "D2P:" << std::endl;
    D2P.print();
    std::cout << "D3P:" << std::endl;
    D3P.print();

    std::cout << "-D1P:" << std::endl;
    D1P.inv().print();

    return 0;
}
