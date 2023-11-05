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

    // 10M, 3S
    Number ZZ = Z1 * Z2;
    Number U11Z2 = U11 * Z2;
    Number U21Z1 = U21 * Z1;
    Number U10Z2 = U10 * Z2;
    Number U20Z1 = U20 * Z1;
    Number V11Z2 = V11 * Z2;
    Number V21Z1 = V21 * Z1;
    Number V10Z2 = V10 * Z2;
    Number V20Z1 = V20 * Z1;

    Number U11Z2S = U11Z2 * U11Z2;
    Number U21Z1S = U21Z1 * U21Z1;

    Number T11 = U11 * U11;

    Number Z1S = Z1 * Z1;
    Number ZZ2 = ZZ * ZZ;

    // 3M

    // (Z1Z2)^2 がかかっている．
    Number M1 = U11Z2S - U21Z1S + ZZ * (U20Z1 - U10Z2);
    Number M2 = U21Z1 * U20Z1 - U11Z2 * U10Z2;
    // Z1Z2 がかかっている．
    Number M3 = U11Z2 - U21Z1;
    Number M4 = U20Z1 - U10Z2;

    Number z1 = V10Z2 - V20Z1;
    Number z2 = V11Z2 - V21Z1;

    // 4M

    // (Z1Z2)^3 がかかっている．
    Number t1 = (M2 - z1) * (z2 - M1);
    Number t2 = (-z1 - M2) * (z2 + M1);
    // (Z1Z2)^2 がかかっている．
    Number t3 = (-z1 + M4) * (z2 - M3);
    Number t4 = (-z1 - M4) * (z2 + M3);

    // 1M
    // (Z1Z2)^3 がかかっている．
    Number l2_num = t1 - t2;
    // (Z1Z2)^2 がかかっている．
    Number l3_num = t3 - t4;
    // (Z1Z2)^3 がかかっている．
    Number d = -(M2 - M4) * (M1 + M3);
    d = d + d + (t3 + t4) - t1 - t2;

    // std::cout << l2_num << " " << l3_num << " " << d << " " << ZZ << std::endl;

    // 2M, 2S
    Number A = d * d;
    Number B = ZZ2 * l3_num * l3_num - f6 * A;

    // 6M
    // (l3_num^2 ZZ^2 - f6 d^2) ZZ^2 がかかっている．
    Number l2l3ZZ = l2_num * l3_num * ZZ;
    Number U1d = - B * (U11Z2 + U21Z1) - (f5 * A - l2l3ZZ - l2l3ZZ) * ZZ;
    Number Zd = B * ZZ;
    
    // 21M, 1S
    Number U0d = l3_num * ZZ2 * (U10 * Z1 - T11) + ZZ * Z1 * (l2_num * U11 + d * V11);
    U0d = (U0d + U0d) * l3_num + Z1S * (l2_num * l2_num - A * f4);
    U0d = Zd * Z2 * U0d;
    U0d = U0d - B * Z1 * ((U10Z2 + U20Z1 + U11 * U21) * Zd + (U11Z2 + U21Z1) * U1d);
    Number ZdM = Z2 * Z1S * B;
    Zd = Zd * ZdM;
    U1d = U1d * ZdM;
    Number ZdS = Zd * Zd;

    // 23M, 1S
    Number dZ1ZdS = d * Z1 * ZdS;
    Number ZZL3 = ZZ * l3_num;
    Number Z1ZdL2 = Z1 * Zd * l2_num;
    Number V1d = -Z1ZdL2 * (U11 * Zd - U1d * Z1) - ZZL3 * ( (U1d * U1d - U0d * Zd) * Z1S - (T11 - U10 * Z1) * ZdS) - V11 * dZ1ZdS;
    Number V0d = Z1ZdL2 * (U0d * Z1 - U10 * Zd) + ZZL3 * (U11 * U10 * ZdS - U1d * U0d * Z1S) - V10 * dZ1ZdS;

    // 5M
    ZdM = Zd * Z1S * d;
    Zd = Zd * ZdM;

    U1d = U1d * ZdM;
    U0d = U0d * ZdM;

    WeightedProjectiveMumford ret(f, h, U1d, U0d, V1d, V0d, Zd1, Zd2);
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
