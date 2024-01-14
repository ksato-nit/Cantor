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

    mpz_init_set_str(ExtendedNumber::CHARA, "84922647928907723895182629923699618983233511623459978930846666887434840654527", 10);
    mpz_init_set_str(ExtendedNumber::MCHARA, "-84922647928907723895182629923699618983233511623459978930846666887434840654527", 10);

    int fc[7] = {-1, 3, 6, -2, -3, 1, 1};
    int hc[1] = {0};

    ExtendedPolynomial f(6, fc);
    ExtendedPolynomial h(0, hc);

    mpz_set_str(f.coeff[6].re, "49486469336020032312363114753237753538045189302709837812889891568789329777106", 10);
    mpz_set_str(f.coeff[5].re, "18146572176684363329466943384936827793858523936395648408139279140336594949430", 10);
    mpz_set_str(f.coeff[4].re, "75520695385453195425619146228325522268537048648879238745586077477241925048911", 10);
    mpz_set_str(f.coeff[3].re, "55661282090622854330881331548845511756521092337654843688125103834062467631270", 10);
    mpz_set_str(f.coeff[2].re, "80833291592081646281545155602592608433782855098617566114028816890992939330947", 10);
    mpz_set_str(f.coeff[1].re, "60698319343619711865647044293933240874784308553876779022618656869472568360653", 10);
    mpz_set_str(f.coeff[0].re, "48384504621689098848513719170178371090837990470273004576471363410971726768086", 10);

    ExtendedNumber u11, u10, v11, v10;
    mpz_set_str(u11.re, "84922647928907723895182629923699618983233511623459978930846666887434840654526", 10);
    mpz_set_str(u11.im, "84922647928907723895182629923699618983233511623459978930846666887434840654525", 10);
    mpz_set_str(u10.re, "0", 10);
    mpz_set_str(u10.im, "2", 10);
    mpz_set_str(v11.re, "-8553659501684724959321568441612056953537199202944754025901731139734305203681", 10);
    mpz_set_str(v11.im, "-78122558991811386149920331427818894831865826201942720814098068851317634963551", 10);
    mpz_set_str(v10.re, "-29239183280718052636051139603574546339187293785166194736479663822642621446825", 10);
    mpz_set_str(v10.im, "-13600177874192675490524596991761448302735370843034516233497196072234411381952", 10);
    ExtendedMumford D1(f, h, u11, u10, v11, v10);
    ExtendedProjectiveMumford D1P(f, h, u11, u10, v11, v10);


    ExtendedNumber num, num2, num3;
    mpz_init_set_str(num.re, "44922647928907723895182629923699618983233511623459978930846666887434840654527", 10);
    mpz_init_set_str(num.im, "54922647928907723895182629923699618983233511623459978930846666887434840654527", 10);
    mpz_init_set_str(num2.re, "64922647928907723895182629923699618983233511623459978930846666887434840654527", 10);
    mpz_init_set_str(num2.im, "74922647928907723895182629923699618983233511623459978930846666887434840654527", 10);

    std::cout << "加算" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 100000; ++i){
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
    for(int i = 0; i < 100000; ++i){
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
    for(int i = 0; i < 100000; ++i){
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
    for(int i = 0; i < 100000; ++i){
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
