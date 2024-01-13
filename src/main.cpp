#include <iostream>
#include <chrono>
#include "number.hpp"
#include "polynomial.hpp"
#include "mumford.hpp"
#include "mumford_projective.hpp"
#include "mumford_weighted_projective.hpp"

int main(){

    /*
    mpz_t test1, test2, test3;
    mpz_init(test1);
    mpz_init_set_si(test2, 8);
    mpz_init_set_si(test3, 31);
    Mumford::constant_invert(test1, test2, test3);
    std::cout << test1 << std::endl;
    return 0;
    */

    std::chrono::system_clock::time_point start, end;

    mpz_init_set_str(Number::CHARA, "16613960161207197506610974848157905611744466278275346794947826509160636299163", 10);
    mpz_init_set_str(Number::MCHARA, "-16613960161207197506610974848157905611744466278275346794947826509160636299163", 10);

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

    // 有限体上の基本演算
    mpz_t num, num2;
    mpz_init_set_str(num, "10613960161207197506610974848157905611744466278275346794947826509160636299163", 10);
    mpz_init_set_str(num2, "8613960161207197506610974848157905611744466278275346794947826509160636299163", 10);

    std::cout << "コンストラクタ" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 10000; ++i){
        mpz_t numc;
        mpz_init2(numc, 2048);
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;

    std::cout << "コピー" << std::endl;
    mpz_t num3;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 10000; ++i){
        mpz_init_set(num3, num);
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;

    std::cout << "加算" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 100000; ++i){
        mpz_add(num, num, num2);
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;
    int time_add = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
    mpz_mod(num, num, Number::CHARA);

    std::cout << "乗算" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 100000; ++i){
        mpz_mul(num3, num, num2);
        mpz_mod(num3, num, Number::CHARA);
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;
    int time_mul = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();

    std::cout << "2乗算" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 100000; ++i){
        mpz_mul(num3, num, num);
        mpz_mod(num3, num, Number::CHARA);
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;
    int time_sqr = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();

    std::cout << "逆元" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 100000; ++i){
        Mumford::constant_invert(num, num2, Number::CHARA);
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;
    int time_inv = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();

    std::cout << "I / M : " << (float) time_inv / time_mul << std::endl;
    std::cout << "M / S : " << (float) time_mul / time_sqr << std::endl;
    std::cout << "M / a : " << (float) time_mul / time_add << std::endl;

    // 超楕円曲線上のスカラー倍算
    mpz_class k(Number::CHARA);
    mpz_class one_class(1);
    k = k + one_class;
    k = k * k;
    //k.set_str("166032509583933056829073411341080509916695444781843190281436838500774011631843644521706667298863907454077624083474972757016101982590047184567811979877299804627575859440070530297919737504772537542984827651623100281300498695268416048516481825853425612212381809485673206382258496996410230458686986871205061771101", 10);

    std::cout << "射影 Lange" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 10; ++i){
        ProjectiveMumford Dk = D1P * k;
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;

    std::cout << "アフィン Lange" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 10; ++i){
        Mumford Dk = D1 * k;
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;

    std::cout << "アフィン Costello" << std::endl;
    start = std::chrono::system_clock::now();
    for(int i = 0; i < 10; ++i){
        Mumford Dk = D1.CostelloScalarMultiple(k);
    }
    end = std::chrono::system_clock::now();
    std::cout << "処理時間:" << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;

    return 0;
}
