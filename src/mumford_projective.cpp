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

ProjectiveMumford::ProjectiveMumford(Polynomial f, Polynomial h, Number U1, Number U0, Number V1, Number V0, Number Z, Number W1, Number W0){
    this->f = f;
    this->h = h;
    this->U1 = U1;
    this->U0 = U0;
    this->V1 = V1;
    this->V0 = V0;
    this->Z = Z;
    this->W1 = W1;
    this->W0 = W0;
}

ProjectiveMumford::ProjectiveMumford(Polynomial f, Polynomial h, Number U1, Number U0, Number V1, Number V0, Number Z){
    this->f = f;
    this->h = h;
    this->U1 = U1;
    this->U0 = U0;
    this->V1 = V1;
    this->V0 = V0;
    this->Z = Z;
    this->W1 = U1 * U1;
    this->W0 = U1 * U0;
}

ProjectiveMumford::ProjectiveMumford(Polynomial f, Polynomial h, Number U1, Number U0, Number V1, Number V0){
    this->f = f;
    this->h = h;
    this->U1 = U1;
    this->U0 = U0;
    this->V1 = V1;
    this->V0 = V0;
    this->Z = Number::ONE();
    this->W1 = U1 * U1;
    this->W0 = U1 * U0;
}

ProjectiveMumford ProjectiveMumford::operator + (const ProjectiveMumford& m) const{
    // deg u1 = deg u2 = 2
    ProjectiveMumford ret = this->CostelloAdd(m);
    return ret;
}

ProjectiveMumford ProjectiveMumford::CostelloAdd(const ProjectiveMumford& m) const{
    std::cout << "Projective Costello Addition." << std::endl;

    // Total: 76M, 9S

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

    ProjectiveMumford ret(f, h, U1d, U0d, V1d, V0d, Zd);
    return ret;
}

ProjectiveMumford ProjectiveMumford::LangeAdd(const ProjectiveMumford& m) const{
    std::cout << "Projective Lange Addition." << std::endl;

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

    // 64M, 6S

    // 1. 終結式を計算．
    // 8M, 2S
    Number z1 = U11 * Z2 - U21 * Z1;
    Number z2 = U20 * Z1 - U10 * Z2;
    Number z3 = U11 * z1 + z2 * Z1;
    Number r = z2 * z3 + z1 * z1 * U10; // Z1^3 Z2^2 がかかっている．
    Number rs = r * r;

    // 2. almost inverse を計算．
    // 6M
    Number inv1 = z1;
    Number inv0 = z3;

    Number w0 = V10 * Z2 - V20 * Z1;
    Number w1 = V11 * Z2 - V21 * Z1;
    Number w2 = inv0 * w0;
    Number w3 = inv1 * w1;

    // 3. s を計算．
    // 4M
    Number s1 = (inv0 + Z1 * inv1) * (w0 + w1) - w2 - w3 * (Z1 + U11);
    Number s0 = w2 - U10 * w3;

    // 4. l を計算．全体に (Z1 Z2)^3 がかかっている．
    // 5M
    Number l3 = s1 * Z2;
    Number l2 = s1 * U21;
    Number l0 = s0 * U20;
    Number l1 = (s1 + s0) * (U21 + U20) - l2 - l0; //s1 * U20 + s0 * U21;
    l2 = l2 + s0 * Z2;

    // 5. U' を計算．全体に (Z1 Z2)^6 がかかっている．
    // 18M, 1S
    Number f5Z2 = f5 * Z2;
    Number f6U21 = f6 * U21;
    Number rV21 = r * V21;

    Number t4 = (s1 * l3 - rs * Z2 * f6) * Z2;
    Number t3 = ((l2 * s1 + l3 * s0) - rs * (f5Z2 - f6U21)) * Z2;
    Number t2 = Z2 * (s0 * l2 + s1 * (l1 + rV21 * 2)) - rs * ( (f4 * Z2 - f6 * U20) * Z2 - (f5Z2 - f6U21) * U21 );
 
    // 8M, 2S
    Number t4U11 = t4 * U11;
    Number t3Z1 = t3 * Z1;
    Number Ud2 = t4;
    Ud2 = Ud2 * Z1 * Z1;
    Number Ud1 = t3Z1 - t4U11;
    Ud1 = Ud1 * Z1;
    Number Ud0 = (t2 * Z1 - t4 * U10) * Z1 - (t3Z1 - t4U11) * U11;
    Number Zd = Ud2;
    Number ZdS = Zd * Zd;

    // 6. V' を計算．
    // 10M, 1S
    Number Vd0 = (-l0 - V20 * r) * ZdS - Ud0 * (Ud1 * l3 - l2 * Zd);
    Number Vd1 = (-l1 - rV21) * ZdS - l3 * (Ud1 * Ud1 - Ud0 * Zd) + Zd * Ud1 * l2;

    // 7. Z' を調整．
    // 5M
    Number M = Zd * Z2 * r;
    Ud1 = Ud1 * M;
    Ud0 = Ud0 * M;
    Zd = Zd * M;

    ProjectiveMumford ret(f, h, Ud1, Ud0, Vd1, Vd0, Zd);
    return ret;
}

ProjectiveMumford ProjectiveMumford::doubling() const{
    std::cout << "Projective Lange Doubling." << std::endl;

    Number U1 = this->U1;
    Number U0 = this->U0;

    Number V1 = this->V1;
    Number V0 = this->V0;

    Number Z = this->Z;

    Number f6 = this->f.coeff[6];
    Number f5 = this->f.coeff[5];
    Number f4 = this->f.coeff[4];
    Number f3 = this->f.coeff[3];
    Number f2 = this->f.coeff[2];

    // 59M, 9S

    // 1. precomputation.
    // 3M, 2S

    Number U0Z = U0 * Z;
    Number f5Z = f5 * Z;
    Number Z2 = Z * Z;
    Number Z3 = Z2 * Z;
    Number Z4 = Z2 * Z2;

    // 2. v~ と u の終結式を計算．
    // 4M, 2S
    Number V1t = V1 * 2;
    Number V0t = V0 * 2;

    Number W0 = V1 * V1;
    Number U1s = U1 * U1;
    Number W1 = U1s;
    Number W2 = W0 * 4;
    Number W3 = U1 * V1t;
    // Z^3 がかかっている．
    Number V0tZ = V0t * Z;
    Number r = U0 * W2 + V0t * (V0tZ - W3);
    std::cout << Z << std::endl;
    std::cout << r << std::endl;

    // 3. almost inverse を計算．
    // Z がかかっている．
    Number inv1d = -V1t;
    // Z^2 がかかっている．
    Number inv0d = V0tZ - W3;

    // 4. k を計算．
    // 15M
    Number k4 = f6;
    Number k4U1 = k4 * U1;
    Number k4U0Z = k4 * U0Z;
    Number k3 = f5Z - k4U1;
    Number k3U0Z = k3 * U0Z;
    Number f4Z2 = f4 * Z2;
    Number k2 = f4Z2 - k4U0Z - k3 * U1;
    Number k1 = f3 * Z3 - k3U0Z - k2 * U1;
    Number k0 = f2 * Z4 - W0 * Z2 - k2 * U0Z - k1 * U1;
    // Z^3 がかかっている．
    Number k1d = k1 + W1 * (k3 - k4U1) - k3U0Z + U1 * (k4U0Z * 2 - k2);
    // Z^4 がかかっている．
    Number k0d = k0 + U0Z * (U1 * (k3 - k4U1) + k4U0Z - k2);
    std::cout << k4 << " " << k3 << " " << k2 << " " << k1 << " " << k0 << " " << k1d << " " << k0d << std::endl;

    // 5. s を計算．
    // 5M
    W0 = k0d * inv0d;
    W1 = k1d * inv1d;
    // Z^5 がかかっている．
    Number s1d = (inv0d + inv1d) * (k0d + k1d) - W0 - W1 * (Number::ONE() + U1);
    // Z^6 がかかっている．
    Number s0d = W0 - W1 * U0Z;
    std::cout << s1d << " " << s0d << std::endl;

    // 6. U' を計算．
    // 12M, 4S
    Number rs = r * r;
    Number f5Zk = f5Z - k4U1 * 2;
    Number Z3r = Z3 * r;
    Number rsZ4 = rs * Z4;
    // Z^10 がかかっている．
    Number Ud2 = s1d * s1d - rsZ4 * f6;
    // Z^11 がかかっている．
    Number Ud1 = s1d * s0d * 2 - rsZ4 * f5Zk;
    // Z^12 がかかっている．
    Number V1Z3r = V1 * Z3r;
    Number Ud0 = s0d * s0d + s1d * V1Z3r * 2 - rsZ4 * (f4Z2 - (U0Z * 2 + U1s) * f6 - U1 * f5Zk * 2);
    Ud1 = Ud1 * Z;
    Number Zd = Ud2 * Z2;
    Number Zd2 = Zd * Zd;

    // 7. V' を計算．
    // 15M, 1S
    // r Z^7 がかかっている．
    Number l3 = s1d * Z2;
    Number l2 = (s1d * U1 + s0d) * Z;
    Number l1 = (s1d * U0Z + s0d * U1);
    Number l0 = s0d * U0;
    Number l2Zd = l2 * Zd;
    std::cout << l3 << " " << l2 << " " << l1 << " " << l0 << std::endl;
    Number Vd1 = l3 * (Ud1 * Ud1 - Ud0 * Zd) - Ud1 * l2Zd + (l1 + V1Z3r) * Zd2;
    Number Vd0 = (l3 * Ud1 - l2Zd) * Ud0 + (l0 + V0 * Z3r) * Zd2;

    // 8. 調整
    // 5M
    Number ZM = Z4 * Zd * r;
    Ud1 = Ud1 * ZM;
    Ud0 = Ud0 * ZM;
    Zd = Zd * ZM;
    std::cout << Vd1 << " " << Vd0 << " " << Zd << std::endl;


    ProjectiveMumford ret(f, h, Ud1, Ud0, Vd1, Vd0, Zd);
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

    Polynomial u(2); Polynomial v(1);
    u.coeff[2] = Number::ONE(); u.coeff[1] = this->U1 / this->Z; u.coeff[0] = this->U0 / this->Z;
    v.coeff[1] = this->V1 / this->Z; v.coeff[0] = this->V0 / this->Z;

    std::cout << "[" << u << ", " << v << "]" << std::endl;
    
    return;
}
