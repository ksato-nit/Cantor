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

    mpz_init_set_str(ExtendedNumber::CHARA, "827512402471234900661632350465842959144524719813594547964188836612988107458951679162314438739350392195712435445958395466158429317340203911888824466573450079820676517044845033418304624062441761062718649659756375120218401112241165371", 10);
    mpz_init_set_str(ExtendedNumber::MCHARA, "-827512402471234900661632350465842959144524719813594547964188836612988107458951679162314438739350392195712435445958395466158429317340203911888824466573450079820676517044845033418304624062441761062718649659756375120218401112241165371", 10);

    int fc[7] = {-1, 3, 6, -2, -3, 1, 1};
    int hc[1] = {0};

    ExtendedPolynomial f(6, fc);
    ExtendedPolynomial h(0, hc);

    mpz_set_str(f.coeff[6].re, "408574927522568260646494190818708344923289394155246733745057061347760234777470744844105873004280234841872260152482842131295379950502171648756403167430258933763756867838745876407324986689616739945861628369646848583019500730413413753", 10);
    mpz_set_str(f.coeff[5].re, "646057604203333654565256199725133766858940477418311967251502430873928928969993054398880539279033711849180581910560923929784422762246976366008190607196936483836481028101473110360898357421566961591191720911698797969718903542235882367", 10);
    mpz_set_str(f.coeff[4].re, "614846262136796668250655182767419943446862193747402912871652070113162719458104150784643841966874312233503425991457853010573970036996279554140411488719755857254798210214453124876936218220833154347644956217046905577687168956534810932", 10);
    mpz_set_str(f.coeff[3].re, "44522322012461199512858089421650491553172038916456638692823617727365618822036792482337860498020502189482770385012328307462081813997242322524299904163848837661698821940466336870590251154454717342689786345725916838070074458110815582", 10);
    mpz_set_str(f.coeff[2].re, "2663399975623008568771749256434924933712190070176531983766239622293977948685545379088464299476705912240779667876964288477991932399534862098611275954725657398136434378965065920455563128373533597035484747846700016450949895167871434", 10);
    mpz_set_str(f.coeff[1].re, "253735127427402269957171415195144193555348677146048048977754098272988896681569702318897457456551312141302247202496442852799188159579340786644094532787119343307335381885493212837741767438120836858991415141168905245555967267061165319", 10);
    mpz_set_str(f.coeff[0].re, "182263745344618790831340607620523635307141339303832513911263369992129325482652613017488724081311772352944481837661799097840061256535133200419434575036392187211782384474967980229921555090108308381413079360433791281989721333292097365", 10);

    ExtendedNumber u11, u10, v11, v10;
    mpz_set_str(u11.re, "827512402471234900661632350465842959144524719813594547964188836612988107458951679162314438739350392195712435445958395466158429317340203911888824466573450079820676517044845033418304624062441761062718649659756375120218401112241165368", 10);
    mpz_set_str(u11.im, "827512402471234900661632350465842959144524719813594547964188836612988107458951679162314438739350392195712435445958395466158429317340203911888824466573450079820676517044845033418304624062441761062718649659756375120218401112241165370", 10);
    mpz_set_str(u10.re, "2", 10);
    mpz_set_str(u10.im, "2", 10);
    mpz_set_str(v11.re, "-688369786586276896261828266196434108056662915631039836284521335287711090466659873513494822113957852181324935504139137192151565149474075762851030669610177953832039288398189561024109202314737554964748552710606105697857575280129751986", 10);
    mpz_set_str(v11.im, "-764464472490077970803528340562694704445436617711900898293215936776715728153253295616071656320040222885753711744740684847389020544288163843961136253307338280189484157520778670693248802718738030175952478430391024130688650882050734289", 10);
    mpz_set_str(v10.re, "-54723679893387964512678971107857693878740270357564635721333418934436513493179630943686886835418643119788794481067987474231700061319853789099925937837319563481927302986079350974615616947481689590558426266160545725686306587354470301", 10);
    mpz_set_str(v10.im, "-126095859962313859716208019806296509398176204203387299341945799672544758611396767092485564838620338619917447402435421237538817546104080135855376426532223599262384719048132725450111642687407461773532342458730701979059500460380862164", 10);
    ExtendedMumford D1(f, h, u11, u10, v11, v10);
    ExtendedProjectiveMumford D1P(f, h, u11, u10, v11, v10);


    ExtendedNumber num, num2, num3;
    mpz_init_set_str(num.re, "327512402471234900661632350465842959144524719813594547964188836612988107458951679162314438739350392195712435445958395466158429317340203911888824466573450079820676517044845033418304624062441761062718649659756375120218401112241165371", 10);
    mpz_init_set_str(num.im, "427512402471234900661632350465842959144524719813594547964188836612988107458951679162314438739350392195712435445958395466158429317340203911888824466573450079820676517044845033418304624062441761062718649659756375120218401112241165371", 10);
    mpz_init_set_str(num2.re, "527512402471234900661632350465842959144524719813594547964188836612988107458951679162314438739350392195712435445958395466158429317340203911888824466573450079820676517044845033418304624062441761062718649659756375120218401112241165371", 10);
    mpz_init_set_str(num2.im, "627512402471234900661632350465842959144524719813594547964188836612988107458951679162314438739350392195712435445958395466158429317340203911888824466573450079820676517044845033418304624062441761062718649659756375120218401112241165371", 10);

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
