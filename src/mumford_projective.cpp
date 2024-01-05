#include "mumford_projective.hpp"

ProjectiveMumford::ProjectiveMumford(){
    this->f = Polynomial();
    this->h = Polynomial();
    this->U1 = Number::ZERO();
    this->U0 = Number::ONE();
    this->V1 = Number::ZERO();
    this->V0 = Number::ZERO();
    this->Z = Number::ONE();
}

ProjectiveMumford::ProjectiveMumford(Polynomial f, Polynomial h){
    this->f = f;
    this->h = h;
    this->U1 = Number::ZERO();
    this->U0 = Number::ONE();
    this->V1 = Number::ZERO();
    this->V0 = Number::ZERO();
    this->Z = Number::ONE();
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

ProjectiveMumford ProjectiveMumford::operator * (const mpz_class& k_) const{
    Polynomial f = this->f;
    Polynomial h = this->h;
    mpz_class k = k_;

    // double-and-add method によりスカラー倍を計算する．

    // 連除法で 2 進数に変換．
    mpz_class two = 2;
    int count = 0;
    std::vector<int> bits;
    while(true){
        mpz_class rem = k % 2;
        bits.push_back((int) rem.get_si());

        k = k / 2;
        ++count;

        if(k < 1){
            break;
        }
    }

    ProjectiveMumford D = ProjectiveMumford::zero(f, h);
    ProjectiveMumford now = *this;
    for(int i = 0; i < count; ++i){
        if(bits[i] == 1){
            if(D.isZero()){
                D = now;
            }else if(now.isZero()){
                // D = D;
            }else{
                D = D.LangeAdd(now);
            }
            //D.print();
        }
        now = now.LangeDoubling();
    }
    return D;
}

ProjectiveMumford ProjectiveMumford::CostelloAdd(const ProjectiveMumford& m) const{
    //std::cout << "Projective Costello Addition." << std::endl;

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
    //std::cout << "Projective Lange Addition." << std::endl;

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

    Number temp1, temp2;
    // 1. 終結式を計算．
    // 8M, 2S
    Number z1, z2, z3, r, rs;
    mpz_mul(z1.value, U11.value, Z2.value);
    mpz_mod(z1.value, z1.value, Number::CHARA);
    mpz_mul(temp1.value, U21.value, Z1.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);
    mpz_sub(z1.value, z1.value, temp1.value);
    mpz_mul(z2.value, U20.value, Z1.value);
    mpz_mod(z2.value, z2.value, Number::CHARA);
    mpz_mul(temp1.value, U10.value, Z2.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);
    mpz_sub(z2.value, z2.value, temp1.value);
    mpz_mul(z3.value, U11.value, z1.value);
    mpz_mod(z3.value, z3.value, Number::CHARA);
    mpz_mul(temp1.value, z2.value, Z1.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);
    mpz_add(z3.value, z3.value, temp1.value);
    mpz_mul(r.value, z1.value, z1.value);
    mpz_mod(r.value, r.value, Number::CHARA);
    mpz_mul(r.value, r.value, U10.value);
    mpz_mod(r.value, r.value, Number::CHARA);
    mpz_mul(temp1.value, z2.value, z3.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);
    mpz_add(r.value, temp1.value, r.value);
    mpz_mul(rs.value, r.value, r.value);
    mpz_mod(rs.value, rs.value, Number::CHARA);

    // 2. almost inverse を計算．
    // 6M
    //Number inv1 = z1;
    //Number inv0 = z3;
    Number w0, w1, w2, w3;
    mpz_mul(w0.value, V10.value, Z2.value);
    mpz_mod(w0.value, w0.value, Number::CHARA);
    mpz_mul(temp1.value, V20.value, Z1.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);
    mpz_sub(w0.value, w0.value, temp1.value);
    mpz_mul(w1.value, V11.value, Z2.value);
    mpz_mod(w1.value, w1.value, Number::CHARA);
    mpz_mul(temp1.value, V21.value, Z1.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);
    mpz_sub(w1.value, w1.value, temp1.value);
    mpz_mul(w2.value, z3.value, w0.value);
    mpz_mod(w2.value, w2.value, Number::CHARA);
    mpz_mul(w3.value, z1.value, w1.value);
    mpz_mod(w3.value, w3.value, Number::CHARA);

    // 3. s を計算．
    // 4M
    Number s1, s0;
    mpz_add(s1.value, w0.value, w1.value);
    mpz_mul(temp1.value, Z1.value, z1.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);
    mpz_add(temp1.value, temp1.value, z3.value);
    mpz_mul(s1.value, s1.value, temp1.value);
    mpz_mod(s1.value, s1.value, Number::CHARA);
    mpz_sub(s1.value, s1.value, w2.value);
    mpz_add(temp1.value, Z1.value, U11.value);
    mpz_mul(temp1.value, w3.value, temp1.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);
    mpz_sub(s1.value, s1.value, temp1.value);
    mpz_mul(s0.value, U10.value, w3.value);
    mpz_mod(s0.value, s0.value, Number::CHARA);
    mpz_sub(s0.value, w2.value, s0.value);

    // 4. l を計算．全体に (Z1 Z2)^3 がかかっている．
    // 5M
    Number l3, l2, l1, l0;
    mpz_mul(l3.value, s1.value, Z2.value);
    mpz_mod(l3.value, l3.value, Number::CHARA);
    mpz_mul(l2.value, s1.value, U21.value);
    mpz_mod(l2.value, l2.value, Number::CHARA);
    mpz_mul(l0.value, s0.value, U20.value);
    mpz_mod(l0.value, l0.value, Number::CHARA);
    mpz_add(l1.value, s1.value, s0.value);
    mpz_add(temp1.value, U21.value, U20.value);
    mpz_mul(l1.value, l1.value, temp1.value);
    mpz_mod(l1.value, l1.value, Number::CHARA);
    mpz_sub(l1.value, l1.value, l2.value);
    mpz_sub(l1.value, l1.value, l0.value);
    mpz_mul(temp1.value, s0.value, Z2.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);
    mpz_add(l2.value, temp1.value, l2.value);

    // 5. U' を計算．全体に (Z1 Z2)^6 がかかっている．
    // 18M, 1S
    Number f5Z2, f6U21, rV21;
    mpz_mul(f5Z2.value, f5.value, Z2.value);
    mpz_mod(f5Z2.value, f5Z2.value, Number::CHARA);
    mpz_mul(f6U21.value, f6.value, U21.value);
    mpz_mod(f6U21.value, f6U21.value, Number::CHARA);
    mpz_mul(rV21.value, r.value, V21.value);
    mpz_mod(rV21.value, rV21.value, Number::CHARA);

    Number t4, t3, t2;
    mpz_mul(t4.value, rs.value, Z2.value);
    mpz_mod(t4.value, t4.value, Number::CHARA);
    mpz_mul(t4.value, t4.value, f6.value);
    mpz_mod(t4.value, t4.value, Number::CHARA);
    mpz_mul(temp1.value, s1.value, l3.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);
    mpz_sub(t4.value, temp1.value, t4.value);
    mpz_mul(t4.value, t4.value, Z2.value);
    mpz_mod(t4.value, t4.value, Number::CHARA);

    mpz_mul(t3.value, l2.value, s1.value);
    mpz_mod(t3.value, t3.value, Number::CHARA);
    mpz_mul(temp1.value, l3.value, s0.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);
    mpz_add(t3.value, t3.value, temp1.value);
    mpz_sub(temp1.value, f5Z2.value, f6U21.value);
    mpz_mul(temp1.value, temp1.value, rs.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);
    mpz_sub(t3.value, t3.value, temp1.value);
    mpz_mul(t3.value, t3.value, Z2.value);
    mpz_mod(t3.value, t3.value, Number::CHARA);

    mpz_add(t2.value, l1.value, rV21.value);
    mpz_add(t2.value, t2.value, rV21.value);
    mpz_mul(t2.value, t2.value, s1.value);
    mpz_mod(t2.value, t2.value, Number::CHARA);
    mpz_mul(temp1.value, s0.value, l2.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);
    mpz_add(t2.value, t2.value, temp1.value);
    mpz_mul(t2.value, t2.value, Z2.value);
    mpz_mod(t2.value, t2.value, Number::CHARA);
    mpz_mul(temp1.value, f4.value, Z2.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);
    mpz_mul(temp2.value, f6.value, U20.value);
    mpz_mod(temp2.value, temp2.value, Number::CHARA);
    mpz_sub(temp1.value, temp1.value, temp2.value);
    mpz_mul(temp1.value, temp1.value, Z2.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);
    mpz_sub(temp2.value, f5Z2.value, f6U21.value);
    mpz_mul(temp2.value, U21.value, temp2.value);
    mpz_mod(temp2.value, temp2.value, Number::CHARA);
    mpz_sub(temp1.value, temp1.value, temp2.value);
    mpz_mul(temp1.value, temp1.value, rs.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);
    mpz_sub(t2.value, t2.value, temp1.value);
 
    // 8M, 2S
    Number t4U11, t3Z1, Ud2, Ud1, Ud0, ZdS;
    mpz_mul(t4U11.value, t4.value, U11.value);
    mpz_mod(t4U11.value, t4U11.value, Number::CHARA);
    mpz_mul(t3Z1.value, t3.value, Z1.value);
    mpz_mod(t3Z1.value, t3Z1.value, Number::CHARA);
    mpz_mul(Ud2.value, t4.value, Z1.value);
    mpz_mod(Ud2.value, Ud2.value, Number::CHARA);
    mpz_mul(Ud2.value, Ud2.value, Z1.value);
    mpz_mod(Ud2.value, Ud2.value, Number::CHARA);
    mpz_sub(Ud1.value, t3Z1.value, t4U11.value);
    mpz_mul(Ud1.value, Ud1.value, Z1.value);
    mpz_mod(Ud1.value, Ud1.value, Number::CHARA);
    mpz_mul(Ud0.value, t2.value, Z1.value);
    mpz_mod(Ud0.value, Ud0.value, Number::CHARA);
    mpz_mul(temp1.value, t4.value, U10.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);
    mpz_sub(Ud0.value, Ud0.value, temp1.value);
    mpz_mul(Ud0.value, Ud0.value, Z1.value);
    mpz_sub(temp1.value, t3Z1.value, t4U11.value);
    mpz_mul(temp1.value, temp1.value, U11.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);
    mpz_sub(Ud0.value, Ud0.value, temp1.value);
    mpz_mul(ZdS.value, Ud2.value, Ud2.value);
    mpz_mod(ZdS.value, ZdS.value, Number::CHARA);

    // 6. V' を計算．
    // 10M, 1S
    Number Vd0, Vd1;
    mpz_mul(Vd0.value, Ud1.value, l3.value);
    mpz_mod(Vd0.value, Vd0.value, Number::CHARA);
    mpz_mul(temp1.value, l2.value, Ud2.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);
    mpz_sub(Vd0.value, Vd0.value, temp1.value);
    mpz_mul(Vd0.value, Vd0.value, Ud0.value);
    mpz_mod(Vd0.value, Vd0.value, Number::CHARA);
    mpz_mul(temp1.value, V20.value, r.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);
    mpz_add(temp1.value, temp1.value, l0.value);
    mpz_mul(temp1.value, temp1.value, ZdS.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);
    mpz_add(Vd0.value, Vd0.value, temp1.value);
    mpz_neg(Vd0.value, Vd0.value);

    mpz_mul(Vd1.value, Ud1.value, Ud1.value);
    mpz_mod(Vd1.value, Vd1.value, Number::CHARA);
    mpz_mul(temp1.value, Ud0.value, Ud2.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);
    mpz_sub(Vd1.value, Vd1.value, temp1.value);
    mpz_mul(Vd1.value, l3.value, Vd1.value);
    mpz_mod(Vd1.value, Vd1.value, Number::CHARA);
    mpz_mul(temp1.value, Ud2.value, Ud1.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);
    mpz_mul(temp1.value, temp1.value, l2.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);
    mpz_sub(Vd1.value, temp1.value, Vd1.value);
    mpz_add(temp1.value, l1.value, rV21.value);
    mpz_mul(temp1.value, temp1.value, ZdS.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);
    mpz_sub(Vd1.value, Vd1.value, temp1.value);

    // 7. Z' を調整．
    // 5M
    Number M, Zd;
    mpz_mul(M.value, Ud2.value, Z2.value);
    mpz_mod(M.value, M.value, Number::CHARA);
    mpz_mul(M.value, M.value, r.value);
    mpz_mod(M.value, M.value, Number::CHARA);
    mpz_mul(Ud1.value, Ud1.value, M.value);
    mpz_mod(Ud1.value, Ud1.value, Number::CHARA);
    mpz_mul(Ud0.value, Ud0.value, M.value);
    mpz_mod(Ud0.value, Ud0.value, Number::CHARA);
    mpz_mul(Zd.value, Ud2.value, M.value);
    mpz_mod(Zd.value, Zd.value, Number::CHARA);

    ProjectiveMumford ret(f, h, Ud1, Ud0, Vd1, Vd0, Zd);
    return ret;
}

ProjectiveMumford ProjectiveMumford::LangeDoubling() const{
    //std::cout << "Projective Lange Doubling." << std::endl;

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
    Number V1t = V1 + V1;
    Number V0t = V0 + V0;

    Number W0 = V1 * V1;
    Number U1s = U1 * U1;
    Number W1 = U1s;
    Number W2 = W0 + W0;
    W2 = W2 + W2;
    Number W3 = U1 * V1t;
    // Z^3 がかかっている．
    Number V0tZ = V0t * Z;
    Number r = U0 * W2 + V0t * (V0tZ - W3);

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
    Number k1d = k1 + W1 * (k3 - k4U1) - k3U0Z + U1 * (k4U0Z + k4U0Z - k2);
    // Z^4 がかかっている．
    Number k0d = k0 + U0Z * (U1 * (k3 - k4U1) + k4U0Z - k2);

    // 5. s を計算．
    // 5M
    W0 = k0d * inv0d;
    W1 = k1d * inv1d;
    // Z^5 がかかっている．
    Number s1d = (inv0d + inv1d) * (k0d + k1d) - W0 - W1 * (Number::ONE() + U1);
    // Z^6 がかかっている．
    Number s0d = W0 - W1 * U0Z;

    // 6. U' を計算．
    // 12M, 4S
    Number rs = r * r;
    Number f5Zk = f5Z - k4U1 - k4U1;
    Number Z3r = Z3 * r;
    Number rsZ4 = rs * Z4;
    // Z^10 がかかっている．
    Number Ud2 = s1d * s1d - rsZ4 * f6;
    // Z^11 がかかっている．
    Number Ud1 = s1d * s0d;
    Ud1 = Ud1 + Ud1;
    Ud1 = Ud1 - rsZ4 * f5Zk;
    // Z^12 がかかっている．
    Number V1Z3r = V1 * Z3r;
    Number Ud0 = s0d * s0d + s1d * V1Z3r * 2 - rsZ4 * (f4Z2 - (U0Z + U0Z + U1s) * f6 - U1 * f5Zk * 2);
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
    Number Vd1 = -l3 * (Ud1 * Ud1 - Ud0 * Zd) + Ud1 * l2Zd - (l1 + V1Z3r) * Zd2;
    Number Vd0 = -(l3 * Ud1 - l2Zd) * Ud0 - (l0 + V0 * Z3r) * Zd2;

    // 8. 調整
    // 5M
    Number ZM = Z4 * Zd * r;
    Ud1 = Ud1 * ZM;
    Ud0 = Ud0 * ZM;
    Zd = Zd * ZM;


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

ProjectiveMumford ProjectiveMumford::zero() const{
    Polynomial f = this->f;
    Polynomial h = this->h;

    ProjectiveMumford zero(f, h);
    return zero;
}

ProjectiveMumford ProjectiveMumford::zero(const Polynomial& f, const Polynomial& h){
    ProjectiveMumford zero(f, h);
    return zero;
}

void ProjectiveMumford::print() const{
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

bool ProjectiveMumford::isZero() const{
    Number U1 = this->U1;
    Number U0 = this->U0;
    Number V1 = this->V1;
    Number V0 = this->V0;
    Number Z = this->Z;

    if(V1.isZero() && V0.isZero()){
        if(U1.isZero() && U0 == Z){
            return true;
        }
    }
    return false;
}
