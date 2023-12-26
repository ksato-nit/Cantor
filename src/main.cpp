#include <iostream>
#include <chrono>
#include "number.hpp"
#include "polynomial.hpp"
#include "mumford.hpp"
#include "mumford_projective.hpp"
#include "mumford_weighted_projective.hpp"

int main(){
    std::chrono::system_clock::time_point start, end;

    Number::CHARA.set_str("16613960161207197506610974848157905611744466278275346794947826509160636299163", 10);
    Number::MCHARA.set_str("-16613960161207197506610974848157905611744466278275346794947826509160636299163", 10);

    int fc[7] = {-1, 3, 6, -2, -3, 1, 1};
    int hc[1] = {0};

    Polynomial f(6, fc);
    Polynomial h(0, hc);

    f.coeff[6].value.set_str("2073000536114376942953123179116934719680183993989471825876059560213989900704", 10);
    f.coeff[5].value.set_str("3530944568397616245686198570421241279588573964842549582651697371642689175556", 10);
    f.coeff[4].value.set_str("10009300617611043175219184591573821675367483763538045027243656482902104726569", 10);
    f.coeff[3].value.set_str("12226457848623687207235732482853906790835753308237257765742579070022899979760", 10);
    f.coeff[2].value.set_str("16328350437044255459278516443823631509653251729098417489738068431685394667846", 10);
    f.coeff[1].value.set_str("1930817213452871091184273173435800486862637979276692503440553590767668851607", 10);
    f.coeff[0].value.set_str("197719815184848422570808706403705284746229851890683993523807438753698171751", 10);

    Number u1, u0, v1, v0;
    u1.value.set_str("16613960161207197506610974848157905611744466278275346794947826509160636299157", 10);
    u0.value = 5;

    v1.value.set_str("-16356435725421449606272205061240827370228669314115845329141561788096218629065", 10);
    v0.value.set_str("-14070429986860860266841156334061046472731180035955518717936741632036392261964", 10);

    Mumford D1(f, h, u1, u0, v1, v0);
    ProjectiveMumford D1P(f, h, u1, u0, v1, v0);

    Number num;
    num.value.set_str("7613960161207197506610974848157905611744466278275346794947826509160636299163", 10);
    Number num2;
    num2.value.set_str("11613960161207197506610974848157905611744466278275346794947826509160636299163", 10);

    // 有限体上の基本演算
    std::cout << "加算" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 1000000; ++i){
        Number num3 = num + num2;
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;

    std::cout << "乗算" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 1000000; ++i){
        Number num3 = num * num2;
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;

    std::cout << "2乗算" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 1000000; ++i){
        Number num3 = num * num;
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;

    std::cout << "逆元" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 1000000; ++i){
        Number num3 = num.inv();
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
