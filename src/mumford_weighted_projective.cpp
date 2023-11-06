#include "mumford_weighted_projective.hpp"

WeightedProjectiveMumford::WeightedProjectiveMumford(){
    this->f = Polynomial();
    this->h = Polynomial();
    this->U1 = Number::ZERO();
    this->U0 = Number::ONE();
    this->V1 = Number::ZERO();
    this->V0 = Number::ZERO();
    this->Z1 = Number::ONE();
    this->Z2 = Number::ONE();
}

WeightedProjectiveMumford::WeightedProjectiveMumford(Polynomial f, Polynomial h){
    this->f = f;
    this->h = h;
    this->U1 = Number::ZERO();
    this->U0 = Number::ONE();
    this->V1 = Number::ZERO();
    this->Z1 = Number::ONE();
    this->Z2 = Number::ONE();
}

WeightedProjectiveMumford::WeightedProjectiveMumford(Polynomial f, Polynomial h, Number U1, Number U0, Number V1, Number V0, Number Z1, Number Z2){
    this->f = f;
    this->h = h;
    this->U1 = U1;
    this->U0 = U0;
    this->V1 = V1;
    this->V0 = V0;
    this->Z1 = Z1;
    this->Z2 = Z2;
}

WeightedProjectiveMumford::WeightedProjectiveMumford(Polynomial f, Polynomial h, Number U1, Number U0, Number V1, Number V0){
    this->f = f;
    this->h = h;
    this->U1 = U1;
    this->U0 = U0;
    this->V1 = V1;
    this->V0 = V0;
    this->Z1 = Number::ONE();
    this->Z2 = Number::ONE();
}

WeightedProjectiveMumford WeightedProjectiveMumford::operator + (const WeightedProjectiveMumford& m) const{
    // deg u1 = deg u2 = 2
    WeightedProjectiveMumford ret = this->CostelloAdd(m);
    return ret;
}

WeightedProjectiveMumford WeightedProjectiveMumford::CostelloAdd(const WeightedProjectiveMumford& m) const{
    std::cout << "Weighted Projective Costello Addition." << std::endl;
    Number U11 = this->U1;
    Number U10 = this->U0;
    Number U21 = m.U1;
    Number U20 = m.U0;

    Number V11 = this->V1;
    Number V10 = this->V0;
    Number V21 = m.V1;
    Number V20 = m.V0;

    Number Z11 = this->Z1;
    Number Z12 = this->Z2;
    Number Z21 = m.Z1;
    Number Z22 = m.Z2;

    Number f6 = this->f.coeff[6];
    Number f5 = this->f.coeff[5];
    Number f4 = this->f.coeff[4];

    Number Z11S = Z11 * Z11;
    Number Z21S = Z21 * Z21;
    Number Z11Q = Z11S * Z11S;
    Number Z21Q = Z21S * Z21S;    
    Number Z12S = Z12 * Z12;
    Number Z22S = Z22 * Z22;
    Number Z11Z21 = Z11 * Z21;
    Number U11Z21S = U11 * Z21S;
    Number U21Z11S = U21 * Z11S;
    Number U10Z21S = U10 * Z21S;
    Number U20Z11S = U20 * Z11S;
    Number V11Z21TZ22 = V11 * Z21S * Z21 * Z22;
    Number V21Z11TZ12 = V21 * Z11S * Z11 * Z12;
    Number V10Z21TZ22 = V10 * Z21S * Z21 * Z22;
    Number V20Z11TZ12 = V20 * Z11S * Z11 * Z12;

    Number ZZAS = Z11S * Z21S;
    Number ZZBS = Z12S * Z22S;

    // (Z11Z21)^4 がかかっている．
    Number M1 = (U11 * U11 - U10 * Z11S) * Z21Q - (U21 * U21 - U20 * Z21S) * Z11Q;
    Number M2 = U21Z11S * U20Z11S - U11Z21S * U10Z21S;
    // (Z11Z21)^3 Z12 Z22 がかかっている．
    Number M3 = U11Z21S - U21Z11S;
    Number M4 = U20Z11S - U10Z21S;

    Number z1 = V10Z21TZ22 - V20Z11TZ12;
    Number z2 = V11Z21TZ22 - V21Z11TZ12;

    // (Z11Z21)^7 Z12 Z22 がかかっている．
    Number t1 = (M2 - z1) * (z2 - M1);
    Number t2 = (-z1 - M2) * (z2 + M1);
    // (Z11Z21)^5 Z12 Z22 がかかっている．
    Number t3 = (-z1 + M4) * (z2 - M3);
    Number t4 = (-z1 - M4) * (z2 + M3);

    // (Z11Z21)^7 Z12 Z22 がかかっている．
    Number l2_num = t1 - t2;
    // (Z11Z21)^5 Z12 Z22 がかかっている．
    Number l3_num = t3 - t4;
    // (Z11Z21)^6 がかかっている．
    Number d = -(M2 - M4) * (M1 + M3);
    d = d + d + (t3 + t4) - t1 - t2;

    Number A = d * d;
    Number B = ZZAS * l3_num * l3_num - f6 * A * ZZBS;

    // (l3_num^2 ZZ^2 - f6 d^2) ZZ^2 がかかっている．
    Number Ud1 = -B * (U11 * Z21S + U21 * Z11S) - (f5 * A * ZZBS - l2_num * l3_num - l2_num * l3_num) * ZZAS;
    Number Zd1 = B * ZZAS;

    std::cout << Ud1 << " " << Zd1 << std::endl;
    
    Number Ud0 = l3_num * (l3_num * Z11S * Z21S * (U10 * Z11S - U11 * U11) * Z12 + (U11 * Z11 * Z12 + V11 * Z11 * Z21 * Z12 * Z22) * Z11) * ZZAS * 2;
    Ud0 = Ud0 + (l2_num * l2_num - f4 * ZZAS * ZZBS * A) * Z11Q * Z12
    Ud0 = Ud0 * Z21 * Zd1;
    Ud0 = Ud0 - (U11 * U21 + (U10 * Z21 * Z21 - U20 * Z11 * Z11) + (U11 * Z21S + U21 * Z11S)) * Z11S * Z12 * Zd1 * (l3_num * Z11Q * Z12Q - f6 * A * ZZAS * Z12 * Z22);

    Number Zd2;
    Number Vd1;
    Number Vd0;

    WeightedProjectiveMumford ret(f, h, Ud1, Ud0, Vd1, Vd0, Zd1, Zd2);
    return ret;
}

WeightedProjectiveMumford WeightedProjectiveMumford::inv(){
    Polynomial f = this->f;
    Polynomial h = this->h;
    Number U1 = this->U1;
    Number U0 = this->U0;
    Number V1 = this->V1;
    Number V0 = this->V0;
    Number Z1 = this->Z1;
    Number Z2 = this->Z2;

    Number h2 = this->h.coeff[2];
    Number h1 = this->h.coeff[1];
    Number h0 = this->h.coeff[0];
    // TODO : ここを書く
    // - (h + v) % u を計算．
    Number V1d;
    Number V0d;
    WeightedProjectiveMumford inv(f, h, U1, U0, V1d, V0d, Z1, Z2);
    return inv;
}

WeightedProjectiveMumford WeightedProjectiveMumford::zero(){
    Polynomial f = this->f;
    Polynomial h = this->h;

    WeightedProjectiveMumford zero(f, h);
    return zero;
}

void WeightedProjectiveMumford::print(){
    std::cout << "[";
    std::cout << this->U1 << " " << this->U0;
    std::cout << ", ";
    std::cout << this->V1 << " " << this->V0;
    std::cout << ", ";
    std::cout << this->Z1 << " " << this->Z2;
    std::cout << "]" << std::endl;
    return;
}
