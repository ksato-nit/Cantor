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
    Mumford sum1 = D1 + D2;
    end = std::chrono::system_clock::now();
    sum1.print();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;

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
    std::cout << "D1WP:" << std::endl;
    D1WP.print();
    std::cout << "D2WP:" << std::endl;
    D2WP.print();

    start = std::chrono::system_clock::now();
    std::cout << "D1WP + D2WP:" << std::endl;
    WeightedProjectiveMumford sum1P = D1WP + D2WP;
    end = std::chrono::system_clock::now();
    sum1P.print();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;

    return 0;
}
