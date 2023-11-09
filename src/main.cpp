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
    int u1c[3] = {5, 25, 1};
    int v1c[2] = {-2, -23};
    int u2c[3] = {17, 17, 1};
    int v2c[2] = {-27, -16};

    Polynomial f(6, fc);
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

    std::chrono::system_clock::time_point start, end;
    
    start = std::chrono::system_clock::now();
    std::cout << "D1 + D2:" << std::endl;
    Mumford sum0 = D1.CantorAdd(D2);
    end = std::chrono::system_clock::now();
    sum0.print();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl << std::endl;

    start = std::chrono::system_clock::now();
    std::cout << "D1 + D2:" << std::endl;
    Mumford sum1 = D1 + D2;
    end = std::chrono::system_clock::now();
    sum1.print();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl << std::endl;

    Number Z1 = Number(2);
    Number Z2 = Number(3);
    Polynomial u1_p = u1 * Z1;
    Polynomial v1_p = v1 * Z1;
    Polynomial u2_p = u2 * Z2;
    Polynomial v2_p = v2 * Z2;

    ProjectiveMumford D1P(f, h, u1_p.coeff[1], u1_p.coeff[0], v1_p.coeff[1], v1_p.coeff[0], Z1);
    ProjectiveMumford D2P(f, h, u2_p.coeff[1], u2_p.coeff[0], v2_p.coeff[1], v2_p.coeff[0], Z2);
    std::cout << "D1P:" << std::endl;
    D1P.print();
    std::cout << "D2P:" << std::endl;
    D2P.print();

    start = std::chrono::system_clock::now();
    std::cout << "D1WP + D2WP:" << std::endl;
    ProjectiveMumford sumP = D1P + D2P;
    end = std::chrono::system_clock::now();
    sumP.print();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl << std::endl;

    Number Z11 = Number(2);
    Number Z12 = Number(3);
    Number Z21 = Number(3);
    Number Z22 = Number(5);
    Polynomial u1_wp = u1 * Z11 * Z11;
    Polynomial v1_wp = v1 * Z11 * Z11 * Z11 * Z12;
    Polynomial u2_wp = u2 * Z21 * Z21;
    Polynomial v2_wp = v2 * Z21 * Z21 * Z21 * Z22;

    WeightedProjectiveMumford D1WP(f, h, u1_wp.coeff[1], u1_wp.coeff[0], v1_wp.coeff[1], v1_wp.coeff[0], Z11, Z12);
    WeightedProjectiveMumford D2WP(f, h, u2_wp.coeff[1], u2_wp.coeff[0], v2_wp.coeff[1], v2_wp.coeff[0], Z21, Z22);
    std::cout << "D1WP:" << std::endl;
    D1WP.print();
    std::cout << "D2WP:" << std::endl;
    D2WP.print();

    start = std::chrono::system_clock::now();
    std::cout << "D1WP + D2WP:" << std::endl;
    WeightedProjectiveMumford sumWP = D1WP + D2WP;
    end = std::chrono::system_clock::now();
    sumWP.print();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl << std::endl;

    return 0;
}
