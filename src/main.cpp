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

    Polynomial f(6, fc);
    Polynomial h(0, hc);

    Polynomial u1(2);
    Polynomial v1(1);
    Polynomial u2(2);
    Polynomial v2(1);

    u1.coeff[2].value.set_str("1", 10);
    u1.coeff[1].value.set_str("9035175640608653763725813625413309605384738499839513553014541853936733378192409004303855749505979464916002548041749121644708560332309999740716480131369041", 10);
    u1.coeff[0].value.set_str("7", 10);

    v1.coeff[1].value.set_str("-7928573549790744462625443879796718405937805359165373178032746950727581671168412107927475798121132844421230118936402354270248459931921328801401227768311323", 10);
    v1.coeff[0].value.set_str("-5973270907783442106295508176286224209946982640978020568613538039586019156959148365294038372717430311738780255106591380721688424640341439734908973776128503", 10);

    u2.coeff[2].value.set_str("1", 10);
    u2.coeff[1].value.set_str("9035175640608653763725813625413309605384738499839513553014541853936733378192409004303855749505979464916002548041749121644708560332309999740716480131369030", 10);
    u2.coeff[0].value.set_str("90", 10);

    v2.coeff[1].value.set_str("-2817892481354397290730295812531016824112991867660765916601721994702454729630238743007991172726760200887385782126312290530216643204811292905718528210918661", 10);
    v2.coeff[0].value.set_str("-6857441259811277883341781773824022848315126553555174862016122998208420207588815269865317852732590391451341674692723790148934805537158312992390250528283055", 10);

    Mumford D1(f, h, u1, v1);
    Mumford D2(f, h, u2, v2);

    std::cout << "D1:" << std::endl;
    D1.print();
    std::cout << "D2:" << std::endl;
    D2.print();

    std::cout << "D1 + D2:" << std::endl;
    Mumford sum1 = D1.LangeAdd(D2);
    sum1.print();

    Polynomial u1_half = u1 * Number(2);
    Polynomial v1_half = v1 * Number(2);
    Polynomial u2_half = u2 * Number(3);
    Polynomial v2_half = v2 * Number(3);

    ProjectiveMumford D1P(f, h, u1_half.coeff[1], u1_half.coeff[0], v1_half.coeff[1], v1_half.coeff[0], Number(2));
    ProjectiveMumford D2P(f, h, u2_half.coeff[1], u2_half.coeff[0], v2_half.coeff[1], v2_half.coeff[0], Number(3));
    std::cout << "D1P:" << std::endl;
    D1P.print();
    std::cout << "D2P:" << std::endl;
    D2P.print();

    std::cout << "D1P + D2P:" << std::endl;
    D1P.LangeAdd(D2P).print();

    return 0;
}
