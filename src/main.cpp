#include <iostream>
#include <chrono>
#include "number.hpp"
#include "polynomial.hpp"
#include "mumford.hpp"
#include "mumford_projective.hpp"
#include "mumford_weighted_projective.hpp"

int main(){
    std::chrono::system_clock::time_point start, end;

    mpz_init_set_str(Number::CHARA, "13291570262456491945649919417472694507895082458159858770373081004488760440120711808381934639404900616533965155060575198692817091309491374216571938755587637", 10);
    mpz_init_set_str(Number::MCHARA, "-13291570262456491945649919417472694507895082458159858770373081004488760440120711808381934639404900616533965155060575198692817091309491374216571938755587637", 10);

    int fc[7] = {-1, 3, 6, -2, -3, 1, 1};
    int hc[1] = {0};

    Polynomial f(6, fc);
    Polynomial h(0, hc);

    mpz_set_str(f.coeff[6].value, "2073000536114376942953123179116934719680183993989471825876059560213989900704", 10);
    mpz_set_str(f.coeff[5].value, "3530944568397616245686198570421241279588573964842549582651697371642689175556", 10);
    mpz_set_str(f.coeff[4].value, "10009300617611043175219184591573821675367483763538045027243656482902104726569", 10);
    mpz_set_str(f.coeff[3].value, "12226457848623687207235732482853906790835753308237257765742579070022899979760", 10);
    mpz_set_str(f.coeff[2].value, "16328350437044255459278516443823631509653251729098417489738068431685394667846", 10);
    mpz_set_str(f.coeff[1].value, "1930817213452871091184273173435800486862637979276692503440553590767668851607", 10);
    mpz_set_str(f.coeff[0].value, "197719815184848422570808706403705284746229851890683993523807438753698171751", 10);

    Number u1, u0, v1, v0;
    mpz_set_str(u1.value, "16613960161207197506610974848157905611744466278275346794947826509160636299157", 10);
    mpz_set_ui(u0.value, 5);

    mpz_set_str(v1.value, "-16356435725421449606272205061240827370228669314115845329141561788096218629065", 10);
    mpz_set_str(v0.value, "-14070429986860860266841156334061046472731180035955518717936741632036392261964", 10);

    Mumford D1(f, h, u1, u0, v1, v0);
    ProjectiveMumford D1P(f, h, u1, u0, v1, v0);

    Number num;
    mpz_set_str(num.value, "7291570262456491945649919417472694507895082458159858770373081004488760440120711808381934639404900616533965155060575198692817091309491374216571938755587637", 10);
    Number num2;
    mpz_set_str(num2.value, "11291570262456491945649919417472694507895082458159858770373081004488760440120711808381934639404900616533965155060575198692817091309491374216571938755587637", 10);
    Number num3;

    mpz_class t1("7291570262456491945649919417472694507895082458159858770373081004488760440120711808381934639404900616533965155060575198692817091309491374216571938755587637");
    mpz_class t2("11291570262456491945649919417472694507895082458159858770373081004488760440120711808381934639404900616533965155060575198692817091309491374216571938755587637");
    mpz_class p("13291570262456491945649919417472694507895082458159858770373081004488760440120711808381934639404900616533965155060575198692817091309491374216571938755587637");
    
    std::cout << "加算" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 1000000; ++i){
        mpz_class t3;
        t3 = t1 + t2;
        t3 = t3 - p;
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;    

    std::cout << "乗算" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 1000000; ++i){
        mpz_class t3;
        t3 = t1 * t2;
        t3 = t3 % p;
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;   

    mpz_t m1, m2;
    mpz_init_set(m1, num.value); mpz_init_set(m2, num2.value);

    // 有限体上の基本演算
    std::cout << "加算" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 1000000; ++i){
        mpz_t m3;
        mpz_init(m3);
        mpz_add(m3, m1, m2);
        if(mpz_cmp(m3, Number::CHARA) > 0){
            mpz_sub(m3, m3, Number::CHARA);
        }
        //num3 = num + num2;
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;

    std::cout << "乗算" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 1000000; ++i){
        mpz_t m3;
        mpz_init(m3);        
        mpz_mul(m3, m1, m2);
        mpz_mod(m3, m3, Number::CHARA);
        //num3 = num * num2;
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;

    std::cout << "2乗算" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 1000000; ++i){
        mpz_t m3;
        mpz_init(m3);        
        mpz_mul(m3, m1, m2);
        mpz_mod(m3, m3, Number::CHARA);
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;

    std::cout << "逆元" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 1000000; ++i){
        mpz_t m3;
        mpz_init(m3);        
        mpz_invert(m3, m2, Number::CHARA);
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;
    return 0;

    // 超楕円曲線上のスカラー倍算
    mpz_class k;
    k.set_str("16613960161207197506610974848157905611744466278275346794947826509160636299164", 10);

    std::cout << "射影 Lange" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 1000; ++i){
        ProjectiveMumford Dk = D1P * k;
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;

    std::cout << "アフィン Lange" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 1000; ++i){
        Mumford Dk = D1 * k;
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;

    std::cout << "アフィン Costello" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 1000; ++i){
        Mumford Dk = D1.CostelloScalarMultiple(k);
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;

    return 0;
}
