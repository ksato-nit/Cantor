#include "mumford_projective.hpp"

ProjectiveMumford::ProjectiveMumford(){
    this->f = Polynomial();
    this->h = Polynomial();
    this->U1 = Number::ZERO();
    this->U0 = Number::ONE();
    this->V1 = Number::ZERO();
    this->V0 = Number::ZERO();
    this->Z = Number::ZERO();
}

ProjectiveMumford::ProjectiveMumford(Polynomial f, Polynomial h){
    this->f = f;
    this->h = h;
    this->U1 = Number::ZERO();
    this->U0 = Number::ONE();
    this->V1 = Number::ZERO();
    this->V0 = Number::ZERO();
    this->Z = Number::ZERO();
}

ProjectiveMumford::ProjectiveMumford(Polynomial f, Polynomial h, Number U1, Number U0, Number V1, Number V0, Number Z){
    this->f = f;
    this->h = h;
    this->U1 = U1;
    this->U0 = U0;
    this->V1 = V1;
    this->V0 = V0;
    this->Z = Z;
}

ProjectiveMumford::ProjectiveMumford(Polynomial f, Polynomial h, Number U1, Number U0, Number V1, Number V0){
    this->f = f;
    this->h = h;
    this->U1 = U1;
    this->U0 = U0;
    this->V1 = V1;
    this->V0 = V0;
    this->Z = Number::ONE();
}

ProjectiveMumford ProjectiveMumford::operator + (const ProjectiveMumford& m) const{
    // deg u1 = deg u2 = 2
    ProjectiveMumford ret = this->CostelloAdd(m);
    return ret;
}

ProjectiveMumford ProjectiveMumford::CostelloAdd(const ProjectiveMumford& m) const{
    std::cout << "Projective Costello Addition." << std::endl;

    Number U11 = this->U1;
    Number U10 = this->U0;
    Number U21 = m.U1;
    Number U20 = m.U0;

    Number V11 = this->V1;
    Number V10 = this->V0;
    Number V21 = m.V1;
    Number V20 = m.V0;

    Number Z1 = this->Z;
    Number Z2 = m.Z;

    Number f6 = this->f.coeff[6];
    Number f5 = this->f.coeff[5];
    Number f4 = this->f.coeff[4];

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
    Number T21 = U21 * U21;

    // (Z1Z2)^2 がかかっている．
    Number M1 = U11Z2S - U21Z1S + ZZ * (U20Z1 - U10Z2);
    Number M2 = U21Z1 * U20Z1 - U11Z2 * U10Z2;
    // Z1Z2 がかかっている．
    Number M3 = U11Z2 - U21Z1;
    Number M4 = U20Z1 - U10Z2;

    std::cout << "M1, M2, M3, M4 :" << std::endl;
    std::cout << M1 << " " << M2 << " " << M3 << " " << M4 << std::endl;

    Number z1 = V10Z2 - V20Z1;
    Number z2 = V11Z2 - V21Z1;

    // (Z1Z2)^4 がかかっている．
    Number t1 = (M2 - ZZ * z1) * (ZZ * z2 - M1);
    Number t2 = (-ZZ * z1 - M2) * (ZZ * z2 + M1);
    // (Z1Z2)^2 がかかっている．
    Number t3 = (-z1 + M4) * (z2 - M3);
    Number t4 = (-z1 - M4) * (z2 + M3);

    // (Z1Z2)^4 がかかっている．
    Number l2_num = t1 - t2;
    // (Z1Z2)^2 がかかっている．
    Number l3_num = t3 - t4;
    // (Z1Z2)^4 がかかっている．
    Number d = -(M2 - ZZ * M4) * (M1 + ZZ * M3) * 2 + (ZZ * ZZ) * (t3 + t4) - t1 - t2;

    std::cout << "t1, t2, t3, t4 :" << std::endl;
    std::cout << t1 << " " << t2 << " " << t3 << " " << t4 << " " << ZZ << " " << ZZ * ZZ << std::endl;
    std::cout << "l3_num, l2_num, d :" << std::endl;
    std::cout << l3_num << " " << l2_num << " " << d << std::endl;

    // ここまで正しい！

    std::cout << "U1, U0 の計算を開始．" << std::endl;

    Number A = d * d;
    Number ZZ2 = ZZ * ZZ;
    Number ZZ4 = ZZ2 * ZZ2;
    Number ZZ6 = ZZ4 * ZZ2;
    Number ZZ8 = ZZ6 * ZZ2;
    Number B = l3_num * l3_num - f6 * A * ZZ4;
    
    /*
    d = d / ZZ4;
    A = d * d;
    Number l3 = l3_num / (ZZ2 * d);
    Number l2 = l2_num / (ZZ4 * d);
    B = l3 * l3 - f6;
    std::cout << l3 << " " << l2 << " " << B << std::endl;
    Number U1d_test = -(B * (U11Z2 + U21Z1) + ZZ * (f5 - l2 * l3 * 2)) / (B * ZZ);
    std::cout << U1d_test << std::endl;
    */

    // (l3_num^2 - f6 d^2 ZZ^4) ZZ^2 がかかっている．
    //std::cout << B << " " << f5 * A * ZZ6 - l2_num * l3_num * 2 << std::endl;
    Number U1d = -B * (U11Z2 + U21Z1) * ZZ - (f5 * A * ZZ6 - l2_num * l3_num * 2);
    Number Zd = B * ZZ2;
    Number Z1S = Z1 * Z1;
    std::cout << U1d << " " << Zd << std::endl;

    Number U0d = l2_num * ZZ * Zd * (l3_num * (U10 * Z1 - T11) + l2_num * U11 * Z1 + d * ZZ * Z1 * V11) * 2;
    U0d = U0d + Z1S * Zd * l2_num - f4 * A * ZZ4 * Z1S * Zd - A * ZZ4 * Z1 * (U1d * Z1 + U10 * Zd) * 2 - A * ZZ4 * Zd * T11;
    Zd = Zd * A * ZZ4 * Z1S;

    ProjectiveMumford ret(f, h);
    return ret;
}

ProjectiveMumford ProjectiveMumford::inv(){
    Polynomial f = this->f;
    Polynomial h = this->h;
    Number U1 = this->U1;
    Number U0 = this->U0;
    Number V1 = this->V1;
    Number V0 = this->V0;
    Number Z = this->Z;

    Number h2 = this->h.coeff[2];
    Number h1 = this->h.coeff[1];
    Number h0 = this->h.coeff[0];
    // - (h + v) % u を計算．
    Number V1d = U1 * h2 - (V1 + h1) * Z;
    Number V0d = U0 * h2 - (V0 + h0) * Z;
    ProjectiveMumford inv(f, h, U1, U0, V1d, V0d, Z);
    return inv;
}

ProjectiveMumford ProjectiveMumford::zero(){
    Polynomial f = this->f;
    Polynomial h = this->h;

    ProjectiveMumford zero(f, h);
    return zero;
}

void ProjectiveMumford::print(){
    std::cout << "[";
    std::cout << this->U1.value << " " << this->U0.value;
    std::cout << ", ";
    std::cout << this->V1.value << " " << this->V0.value;
    std::cout << ", ";
    std::cout << this->Z.value;
    std::cout << "]" << std::endl;
    return;
}
