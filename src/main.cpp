#include <iostream>
#include <chrono>
#include "number.hpp"
#include "extended_number.hpp"
#include "polynomial.hpp"
#include "extended_polynomial.hpp"
#include "mumford.hpp"
#include "extended_mumford.hpp"
#include "extended_mumford_projective.hpp"
#include "mumford_projective.hpp"
#include "mumford_weighted_projective.hpp"

int main(){
    std::chrono::system_clock::time_point start, end;

    mpz_init_set_str(ExtendedNumber::CHARA, "16613960161207197506610974848157905611744466278275346794947826509160636299163", 10);
    mpz_init_set_str(ExtendedNumber::MCHARA, "-16613960161207197506610974848157905611744466278275346794947826509160636299163", 10);

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

    ExtendedNumber u11, u10, v11, v10;
    mpz_set_str(u11.re, "16613960161207197506610974848157905611744466278275346794947826509160636299135", 10);
    mpz_set_str(u11.im, "16613960161207197506610974848157905611744466278275346794947826509160636299163", 10);
    mpz_set_str(u10.re, "52", 10);
    mpz_set_str(u10.im, "18", 10);
    mpz_set_str(v11.re, "-7131143130494287792231338423199247260243653469965222023936866154524651962123", 10);
    mpz_set_str(v11.im, "-8258850595797402173898123732097832688629432892259916631872975170847010583237", 10);
    mpz_set_str(v10.re, "-14986054899981792885389082573241576053098377774554240442752096698406578248938", 10);
    mpz_set_str(v10.im, "-96258969612393158814727383962240234485600493755513531201876167466615132689", 10);
    ExtendedMumford D1(f, h, u11, u10, v11, v10);

    ExtendedProjectiveMumford D1P(f, h, u11, u10, v11, v10);

    mpz_class k;
    k.set_str("16613960161207197506610974848157905611744466278275346794947826509160636299164", 10);

    std::cout << "アフィン Lange" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 1000; ++i){
        ExtendedMumford Dk = D1 * k;
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;

    std::cout << "アフィン Costello" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 1000; ++i){
        ExtendedMumford Dk = D1.CostelloScalarMultiple(k);
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;

    std::cout << "射影Lange" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 1000; ++i){
        ExtendedProjectiveMumford DP = D1P * k;
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;

    return 0;
}
