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

    mpz_init_set_str(ExtendedNumber::CHARA, "20959118493474869068915362541459900434818583513838202077951048445571608878732485727087744553689143906873848018350579", 10);

    int fc[7] = {-1, 3, 6, -2, -3, 1, 1};
    int hc[1] = {0};

    ExtendedPolynomial f(6, fc);
    ExtendedPolynomial h(0, hc);

    mpz_set_str(f.coeff[6].re, "288286611262021659922794675284151922243904163072986219598909498098168808690603954798646070293541684859186025880515", 10);
    mpz_set_str(f.coeff[5].re, "20086690059958422165592436907613121392643418043922104899268771384625963999058874734948178100247888413942604771535144", 10);
    mpz_set_str(f.coeff[4].re, "12809536194993253585749401208004781497306691790580697895256662312542400786519315621721157406153350338712552182872767", 10);
    mpz_set_str(f.coeff[3].re, "10230258218652180276389319373133542533167805569763465127222004142791736142566572581619123209023620589270363777915485", 10);
    mpz_set_str(f.coeff[2].re, "14393897664493686681396987219599161488198538177221187446120927154443042646188208662955081195940829417601495052937604", 10);
    mpz_set_str(f.coeff[1].re, "8030524038373427901574121943483916718524758755254305433414971895021448231073235416059199602300008818523860929656269", 10);
    mpz_set_str(f.coeff[0].re, "20079811380167707324730647988540991996454973616080865582230594061164389570896035371867433850105590218401382181934417", 10);

    ExtendedNumber u11, u10, v11, v10;
    mpz_set_str(u11.re, "20959118493474869068915362541459900434818583513838202077951048445571608878732485727087744553689143906873848018350576", 10);
    mpz_set_str(u11.im, "20959118493474869068915362541459900434818583513838202077951048445571608878732485727087744553689143906873848018350578", 10);
    mpz_set_str(u10.re, "0", 10);
    mpz_set_str(u10.im, "3", 10);
    mpz_set_str(v11.re, "-19613838960309123525435551543410364508295709987746370901835204462790302371616496780885261338951564524105715302330117", 10);
    mpz_set_str(v11.im, "-8987288558586544137700184794615819627029420165043666068723110821375689174833321612423719019632716715652662961357946", 10);
    mpz_set_str(v10.re, "-19109011940007586019687133397080140862676266541981245692539483257812735760037670204507701469111170679601498464499634", 10);
    mpz_set_str(v10.im, "-14956371311190105724730170699072341988548906532545405949732764427016150232965006616904332048480137666789707152627320", 10);
    ExtendedMumford D1(f, h, u11, u10, v11, v10);
    ExtendedProjectiveMumford D1P(f, h, u11, u10, v11, v10);


    ExtendedNumber num, num2, num3;
    mpz_init_set_str(num.re, "44922647928907723895182629923699618983233511623459978930846666887434840654527", 10);
    mpz_init_set_str(num.im, "54922647928907723895182629923699618983233511623459978930846666887434840654527", 10);
    mpz_init_set_str(num2.re, "64922647928907723895182629923699618983233511623459978930846666887434840654527", 10);
    mpz_init_set_str(num2.im, "74922647928907723895182629923699618983233511623459978930846666887434840654527", 10);

    /*
    std::cout << "加算" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 1000000; ++i){
        mpz_add(num.re, num.re, num2.re);
        mpz_add(num.im, num.im, num2.im);
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;
    int time_add = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
    mpz_mod(num.re, num.re, ExtendedNumber::CHARA);
    mpz_mod(num.im, num.im, ExtendedNumber::CHARA);

    mpz_t temp, tempd1, tempd2, tempd3, tempd4;
    mpz_init(temp);
    mpz_init(tempd1);
    mpz_init(tempd2);
    mpz_init(tempd3);
    mpz_init(tempd4);
    ExtendedNumber tempd;
    std::cout << "乗算" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 1000000; ++i){
        mpz_mul(tempd1, num.re, num2.re);
        mpz_mod(tempd1, tempd1, ExtendedNumber::CHARA);
        mpz_mul(tempd2, num.im, num2.im);
        mpz_mod(tempd2, tempd2, ExtendedNumber::CHARA);
        mpz_sub(tempd.re, tempd1, tempd2);
        mpz_add(tempd3, num.re, num.im);
        mpz_add(tempd4, num2.re, num2.im);
        mpz_mul(tempd.im, tempd3, tempd4);
        mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
        mpz_sub(tempd.im, tempd.im, tempd1);
        mpz_sub(tempd.im, tempd.im, tempd2);
        mpz_set(num.re, tempd.re);
        mpz_set(num.im, tempd.im);
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;
    int time_mul = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();

    std::cout << "2乗算" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 1000000; ++i){
        mpz_mul(tempd1, num.re, num.re);
        mpz_mod(tempd1, tempd1, ExtendedNumber::CHARA);
        mpz_mul(tempd2, num.im, num.im);
        mpz_mod(tempd2, tempd2, ExtendedNumber::CHARA);
        mpz_sub(tempd.re, tempd1, tempd2);
        mpz_add(tempd3, num.re, num.im);
        mpz_add(tempd4, num.re, num.im);
        mpz_mul(tempd.im, tempd3, tempd4);
        mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
        mpz_sub(tempd.im, tempd.im, tempd1);
        mpz_sub(tempd.im, tempd.im, tempd2);
        mpz_set(num.re, tempd.re);
        mpz_set(num.im, tempd.im);
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;
    int time_sqr = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();

    std::cout << "逆元" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 1000000; ++i){
        mpz_t denom;
        mpz_init(denom);
        mpz_mul(temp, num.re, num.re);
        mpz_mod(temp, temp, ExtendedNumber::CHARA);
        mpz_mul(denom, num.im, num.im);
        mpz_mod(denom, denom, ExtendedNumber::CHARA);
        mpz_add(denom, denom, temp);
        ExtendedMumford::constant_invert(denom, denom, ExtendedNumber::CHARA);
        mpz_mul(num.re, num.re, denom);
        mpz_mod(num.re, num.re, ExtendedNumber::CHARA);
        mpz_mul(num.im, num.im, denom);
        mpz_mod(num.im, num.im, ExtendedNumber::CHARA);
        mpz_neg(num.im, num.im);
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;
    int time_inv = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();

    std::cout << "I / M : " << (float) time_inv / time_mul << std::endl;
    std::cout << "M / S : " << (float) time_mul / time_sqr << std::endl;
    std::cout << "M / a : " << (float) time_mul / time_add << std::endl;
    */

    mpz_class k(ExtendedNumber::CHARA);
    mpz_class one_class(1);
    k = k + one_class;
    k = k * k;

    std::cout << "アフィン Lange" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 10; ++i){
        ExtendedMumford Dk = D1 * k;
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;

    std::cout << "アフィン Costello" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 10; ++i){
        ExtendedMumford Dk = D1.CostelloScalarMultiple(k);
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;

    std::cout << "射影Lange" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 10; ++i){
        ExtendedProjectiveMumford DP = D1P * k;
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;

    return 0;
}
