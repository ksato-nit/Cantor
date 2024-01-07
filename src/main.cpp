#include <iostream>
#include <chrono>
#include "number.hpp"
#include "extended_number.hpp"
#include "polynomial.hpp"
#include "extended_polynomial.hpp"
#include "mumford.hpp"
#include "extended_mumford.hpp"
#include "mumford_projective.hpp"
#include "mumford_weighted_projective.hpp"

int main(){
    std::chrono::system_clock::time_point start, end;

    mpz_init_set_str(ExtendedNumber::CHARA, "31", 10);
    mpz_init_set_str(ExtendedNumber::MCHARA, "-31", 10);

    int fc[7] = {-1, 3, 6, -2, -3, 1, 1};
    int hc[1] = {0};

    ExtendedPolynomial f(6, fc);
    ExtendedPolynomial h(0, hc);

    mpz_set_str(f.coeff[6].re, "1", 10);
    mpz_set_str(f.coeff[5].re, "1", 10);
    mpz_set_str(f.coeff[4].re, "-3", 10);
    mpz_set_str(f.coeff[3].re, "-2", 10);
    mpz_set_str(f.coeff[2].re, "6", 10);
    mpz_set_str(f.coeff[1].re, "3", 10);
    mpz_set_str(f.coeff[0].re, "-1", 10);

    ExtendedNumber u1, u0, v1, v0;
    mpz_set_str(u1.re, "3", 10);
    mpz_set_str(u1.im, "22", 10);
    mpz_set_str(u0.re, "21", 10);
    mpz_set_str(u0.im, "18", 10);
    mpz_set_str(v1.re, "-3", 10);
    mpz_set_str(v1.im, "-12", 10);
    mpz_set_str(v0.re, "-25", 10);
    mpz_set_str(v0.im, "-8", 10);
    ExtendedMumford D1(f, h, u1, u0, v1, v0);

    D1.print();

    // 超楕円曲線上のスカラー倍算
    ExtendedMumford D = D1.LangeDoubling();
    D.print();

    /*
    std::cout << "アフィン Lange" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 1; ++i){
        ExtendedMumford Dk = D1 * k;
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;

    std::cout << "アフィン Costello" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 1; ++i){
        ExtendedMumford Dk = D1.CostelloScalarMultiple(k);
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;
    */

    return 0;
}
