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

    mpz_init_set_str(ExtendedNumber::CHARA, "166032509583933056829073411341080509916695444781843190281436838500774011631843644521706667298863907454077624083474972757016101982590047184567811979877299804627575859440070530297919737504772537542984827651623100281300498695268416048516481825853425612212381809485673206382258496996410230458686986871205061771101", 10);
    mpz_init_set_str(ExtendedNumber::MCHARA, "-166032509583933056829073411341080509916695444781843190281436838500774011631843644521706667298863907454077624083474972757016101982590047184567811979877299804627575859440070530297919737504772537542984827651623100281300498695268416048516481825853425612212381809485673206382258496996410230458686986871205061771101", 10);

    int fc[7] = {-1, 3, 6, -2, -3, 1, 1};
    int hc[1] = {0};

    ExtendedPolynomial f(6, fc);
    ExtendedPolynomial h(0, hc);

    mpz_set_str(f.coeff[6].re, "52026859858251344126553585298733061327252879662930647106011955807855797347893067797286014619601670129227960270367473227742016371804570180317327881677139767569333280302908148188761701321295454966125292502620575993437234255359864746217794085108089825415376284760328287352436876209059976114600567363828442857084", 10);
    mpz_set_str(f.coeff[5].re, "114140681419114367623251388510658746700913865318451922788166134592709665305268121552599548332803952635310562311157031480628889330477552594572210777686403217185535451919091636755742680745684503112434513528270780602353106315798138078487077277742342339076631631787693058387641891486977183454252020610101528780720", 10);
    mpz_set_str(f.coeff[4].re, "134465360636011097603624967456417316842839570067276874338418859206220699732018959442365530463690422628851113505102991414569825975979700997270531806229485630714123937535309600239542395090019344248216374054296471038574986573657424210834120376311277962539494406130555038648371600869439656245259778190079997770778", 10);
    mpz_set_str(f.coeff[3].re, "58736892402430468027483579871425816145590872434379611382380715023791014601470185172224332119553348875418107789899603794127917000481385778011352678292801847857995307129622710405553287324748249995027573023934071625303715904490097324116677796089875963668398086286958080339501075113692553321238244268150677861723", 10);
    mpz_set_str(f.coeff[2].re, "135381973384021764589641498592103810595048014003487965029082567968423785254583248939933786006269322571787549579471337733402173906084500990375783462813378208416926726823685924645471775767296308718459517477565072076089753414418927494984975981317147865262598682305791521217371555440805578947910733531094483713662", 10);
    mpz_set_str(f.coeff[1].re, "69041738899115788125712924743554391954889901088946525462599079591238351604601677024046981982708239492816763240873184502883972167090732022371345053993158974367040059050661906058935261688526383920272491241805263581578240871657447708093293397206478876807311701068008382676502440713152696356946516470568530409855", 10);
    mpz_set_str(f.coeff[0].re, "90474820555546420017634506750392797998969714290134263938507235429805886467572432252126932526049463493703946746441603563654904088046595650666999442598325048464082792439679061509091378141307715861251065952641577881887805318267537685232386127105149856453681014480720611384214134534293728531236198645166539571099", 10);

    ExtendedNumber u1, u0, v1, v0;
    mpz_set_str(u1.re, "166032509583933056829073411341080509916695444781843190281436838500774011631843644521706667298863907454077624083474972757016101982590047184567811979877299804627575859440070530297919737504772537542984827651623100281300498695268416048516481825853425612212381809485673206382258496996410230458686986871205061771095", 10);
    mpz_set_ui(u0.re, 8);

    mpz_set_str(v1.re, "-93889174465060017039411226397111794564776774281590212285603603485109088963537422689944574440983358254112143711182579114844612669935378076500898468538299402737115517326344562563933249521538114869551103825591588306241858460029064977207840010979742594990081175073907035580312612990212784065944526975182984878819", 10);
    mpz_set_str(v0.re, "-90832142189882902463219062276458138414939494711003979983552055369426342049806688293702707619426167602882924042680513246191647321743107039307876388087617085827842945878705193086861276358832677401263290945158229235690368511239919084994184614037672731320279169041580624690846007235456264204135618681266968241955", 10);

    ExtendedMumford D1(f, h, u1, u0, v1, v0);

    // 超楕円曲線上のスカラー倍算
    mpz_class k;
    k.set_str("166032509583933056829073411341080509916695444781843190281436838500774011631843644521706667298863907454077624083474972757016101982590047184567811979877299804627575859440070530297919737504772537542984827651623100281300498695268416048516481825853425612212381809485673206382258496996410230458686986871205061771101", 10);

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

    return 0;
}
