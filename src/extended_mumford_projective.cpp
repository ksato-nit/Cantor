#include "extended_mumford_projective.hpp"

ExtendedProjectiveMumford::ExtendedProjectiveMumford(){
    this->f = ExtendedPolynomial();
    this->h = ExtendedPolynomial();
    this->U1 = ExtendedNumber::ZERO();
    this->U0 = ExtendedNumber::ONE();
    this->V1 = ExtendedNumber::ZERO();
    this->V0 = ExtendedNumber::ZERO();
    this->Z = ExtendedNumber::ONE();
}

ExtendedProjectiveMumford::ExtendedProjectiveMumford(ExtendedPolynomial f, ExtendedPolynomial h){
    this->f = f;
    this->h = h;
    this->U1 = ExtendedNumber::ZERO();
    this->U0 = ExtendedNumber::ONE();
    this->V1 = ExtendedNumber::ZERO();
    this->V0 = ExtendedNumber::ZERO();
    this->Z = ExtendedNumber::ONE();
}

ExtendedProjectiveMumford::ExtendedProjectiveMumford(ExtendedPolynomial f, ExtendedPolynomial h, ExtendedNumber U1, ExtendedNumber U0, ExtendedNumber V1, ExtendedNumber V0, ExtendedNumber Z, ExtendedNumber W1, ExtendedNumber W0){
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

ExtendedProjectiveMumford::ExtendedProjectiveMumford(ExtendedPolynomial f, ExtendedPolynomial h, ExtendedNumber U1, ExtendedNumber U0, ExtendedNumber V1, ExtendedNumber V0, ExtendedNumber Z){
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

ExtendedProjectiveMumford::ExtendedProjectiveMumford(ExtendedPolynomial f, ExtendedPolynomial h, ExtendedNumber U1, ExtendedNumber U0, ExtendedNumber V1, ExtendedNumber V0){
    this->f = f;
    this->h = h;
    this->U1 = U1;
    this->U0 = U0;
    this->V1 = V1;
    this->V0 = V0;
    this->Z = ExtendedNumber::ONE();
    this->W1 = U1 * U1;
    this->W0 = U1 * U0;
}

ExtendedProjectiveMumford ExtendedProjectiveMumford::operator * (const mpz_class& k_) const{
    ExtendedPolynomial f = this->f;
    ExtendedPolynomial h = this->h;
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

    ExtendedProjectiveMumford D = ExtendedProjectiveMumford::zero(f, h);
    ExtendedProjectiveMumford now = *this;
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

ExtendedProjectiveMumford ExtendedProjectiveMumford::LangeAdd(const ExtendedProjectiveMumford& m) const{
    //std::cout << "Projective Lange Addition." << std::endl;
    // 64M, 6S

    mpz_t temp;
    mpz_init(temp);
    ExtendedNumber temp1, temp2, tempd;
    // 1. 終結式を計算．
    // 8M, 2S
    ExtendedNumber z1, z2, z3, r, rs;
    mpz_mul(z1.re, this->U1.re, m.Z.re);
    mpz_mod(z1.re, z1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->U1.im, m.Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(z1.re, z1.re, temp);
    mpz_mul(z1.im, this->U1.re, m.Z.im);
    mpz_mod(z1.im, z1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->U1.im, m.Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(z1.im, z1.im, temp);



    mpz_mul(temp1.re, m.U1.re, this->Z.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, m.U1.im, this->Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, m.U1.re, this->Z.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, m.U1.im, this->Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_sub(z1.re, z1.re, temp1.re);
    mpz_sub(z1.im, z1.im, temp1.im);
    mpz_mul(z2.re, m.U0.re, this->Z.re);
    mpz_mod(z2.re, z2.re, ExtendedNumber::CHARA);
    mpz_mul(temp, m.U0.im, this->Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(z2.re, z2.re, temp);
    mpz_mul(z2.im, m.U0.re, this->Z.im);
    mpz_mod(z2.im, z2.im, ExtendedNumber::CHARA);
    mpz_mul(temp, m.U0.im, this->Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(z2.im, z2.im, temp);



    mpz_mul(temp1.re, this->U0.re, m.Z.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->U0.im, m.Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, this->U0.re, m.Z.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->U0.im, m.Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_sub(z2.re, z2.re, temp1.re);
    mpz_sub(z2.im, z2.im, temp1.im);
    mpz_mul(z3.re, this->U1.re, z1.re);
    mpz_mod(z3.re, z3.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->U1.im, z1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(z3.re, z3.re, temp);
    mpz_mul(z3.im, this->U1.re, z1.im);
    mpz_mod(z3.im, z3.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->U1.im, z1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(z3.im, z3.im, temp);



    mpz_mul(temp1.re, z2.re, this->Z.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, z2.im, this->Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, z2.re, this->Z.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, z2.im, this->Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_add(z3.re, z3.re, temp1.re);
    mpz_add(z3.im, z3.im, temp1.im);
    mpz_mul(r.re, z1.re, z1.re);
    mpz_mod(r.re, r.re, ExtendedNumber::CHARA);
    mpz_mul(temp, z1.im, z1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(r.re, r.re, temp);
    mpz_mul(r.im, z1.re, z1.im);
    mpz_mod(r.im, r.im, ExtendedNumber::CHARA);
    mpz_mul(temp, z1.im, z1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(r.im, r.im, temp);



    mpz_mul(tempd.re, r.re, this->U0.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, r.im, this->U0.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, r.re, this->U0.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, r.im, this->U0.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(r.re, tempd.re);
    mpz_set(r.im, tempd.im);

    mpz_mul(temp1.re, z2.re, z3.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, z2.im, z3.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, z2.re, z3.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, z2.im, z3.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_add(r.re, temp1.re, r.re);
    mpz_add(r.im, temp1.im, r.im);
    mpz_mul(rs.re, r.re, r.re);
    mpz_mod(rs.re, rs.re, ExtendedNumber::CHARA);
    mpz_mul(temp, r.im, r.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(rs.re, rs.re, temp);
    mpz_mul(rs.im, r.re, r.im);
    mpz_mod(rs.im, rs.im, ExtendedNumber::CHARA);
    mpz_mul(temp, r.im, r.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(rs.im, rs.im, temp);




    // 2. almost inverse を計算．
    // 6M
    //ExtendedNumber inv1 = z1;
    //ExtendedNumber inv0 = z3;
    ExtendedNumber w0, w1, w2, w3;
    mpz_mul(w0.re, this->V0.re, m.Z.re);
    mpz_mod(w0.re, w0.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->V0.im, m.Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(w0.re, w0.re, temp);
    mpz_mul(w0.im, this->V0.re, m.Z.im);
    mpz_mod(w0.im, w0.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->V0.im, m.Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(w0.im, w0.im, temp);



    mpz_mul(temp1.re, m.V0.re, this->Z.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, m.V0.im, this->Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, m.V0.re, this->Z.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, m.V0.im, this->Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_sub(w0.re, w0.re, temp1.re);
    mpz_sub(w0.im, w0.im, temp1.im);
    mpz_mul(w1.re, this->V1.re, m.Z.re);
    mpz_mod(w1.re, w1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->V1.im, m.Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(w1.re, w1.re, temp);
    mpz_mul(w1.im, this->V1.re, m.Z.im);
    mpz_mod(w1.im, w1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->V1.im, m.Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(w1.im, w1.im, temp);



    mpz_mul(temp1.re, m.V1.re, this->Z.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, m.V1.im, this->Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, m.V1.re, this->Z.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, m.V1.im, this->Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_sub(w1.re, w1.re, temp1.re);
    mpz_sub(w1.im, w1.im, temp1.im);
    mpz_mul(w2.re, z3.re, w0.re);
    mpz_mod(w2.re, w2.re, ExtendedNumber::CHARA);
    mpz_mul(temp, z3.im, w0.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(w2.re, w2.re, temp);
    mpz_mul(w2.im, z3.re, w0.im);
    mpz_mod(w2.im, w2.im, ExtendedNumber::CHARA);
    mpz_mul(temp, z3.im, w0.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(w2.im, w2.im, temp);



    mpz_mul(w3.re, z1.re, w1.re);
    mpz_mod(w3.re, w3.re, ExtendedNumber::CHARA);
    mpz_mul(temp, z1.im, w1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(w3.re, w3.re, temp);
    mpz_mul(w3.im, z1.re, w1.im);
    mpz_mod(w3.im, w3.im, ExtendedNumber::CHARA);
    mpz_mul(temp, z1.im, w1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(w3.im, w3.im, temp);




    // 3. s を計算．
    // 4M
    ExtendedNumber s1, s0;
    mpz_add(s1.re, w0.re, w1.re);
    mpz_add(s1.im, w0.im, w1.im);
    mpz_mul(temp1.re, this->Z.re, z1.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->Z.im, z1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, this->Z.re, z1.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->Z.im, z1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_add(temp1.re, temp1.re, z3.re);
    mpz_add(temp1.im, temp1.im, z3.im);
    mpz_mul(tempd.re, s1.re, temp1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s1.im, temp1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, s1.re, temp1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s1.im, temp1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(s1.re, tempd.re);
    mpz_set(s1.im, tempd.im);

    mpz_sub(s1.re, s1.re, w2.re);
    mpz_sub(s1.im, s1.im, w2.im);
    mpz_add(temp1.re, this->Z.re, this->U1.re);
    mpz_add(temp1.im, this->Z.im, this->U1.im);
    mpz_mul(tempd.re, w3.re, temp1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, w3.im, temp1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, w3.re, temp1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, w3.im, temp1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(temp1.re, tempd.re);
    mpz_set(temp1.im, tempd.im);

    mpz_sub(s1.re, s1.re, temp1.re);
    mpz_sub(s1.im, s1.im, temp1.im);
    mpz_mul(s0.re, this->U0.re, w3.re);
    mpz_mod(s0.re, s0.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->U0.im, w3.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(s0.re, s0.re, temp);
    mpz_mul(s0.im, this->U0.re, w3.im);
    mpz_mod(s0.im, s0.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->U0.im, w3.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(s0.im, s0.im, temp);



    mpz_sub(s0.re, w2.re, s0.re);
    mpz_sub(s0.im, w2.im, s0.im);

    // 4. l を計算．全体に (Z1 Z2)^3 がかかっている．
    // 5M
    ExtendedNumber l3, l2, l1, l0;
    mpz_mul(l3.re, s1.re, m.Z.re);
    mpz_mod(l3.re, l3.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s1.im, m.Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(l3.re, l3.re, temp);
    mpz_mul(l3.im, s1.re, m.Z.im);
    mpz_mod(l3.im, l3.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s1.im, m.Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(l3.im, l3.im, temp);



    mpz_mul(l2.re, s1.re, m.U1.re);
    mpz_mod(l2.re, l2.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s1.im, m.U1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(l2.re, l2.re, temp);
    mpz_mul(l2.im, s1.re, m.U1.im);
    mpz_mod(l2.im, l2.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s1.im, m.U1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(l2.im, l2.im, temp);



    mpz_mul(l0.re, s0.re, m.U0.re);
    mpz_mod(l0.re, l0.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s0.im, m.U0.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(l0.re, l0.re, temp);
    mpz_mul(l0.im, s0.re, m.U0.im);
    mpz_mod(l0.im, l0.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s0.im, m.U0.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(l0.im, l0.im, temp);



    mpz_add(l1.re, s1.re, s0.re);
    mpz_add(l1.im, s1.im, s0.im);
    mpz_add(temp1.re, m.U1.re, m.U0.re);
    mpz_add(temp1.im, m.U1.im, m.U0.im);
    mpz_mul(tempd.re, l1.re, temp1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l1.im, temp1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, l1.re, temp1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l1.im, temp1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(l1.re, tempd.re);
    mpz_set(l1.im, tempd.im);

    mpz_sub(l1.re, l1.re, l2.re);
    mpz_sub(l1.im, l1.im, l2.im);
    mpz_sub(l1.re, l1.re, l0.re);
    mpz_sub(l1.im, l1.im, l0.im);
    mpz_mul(temp1.re, s0.re, m.Z.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s0.im, m.Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, s0.re, m.Z.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s0.im, m.Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_add(l2.re, temp1.re, l2.re);
    mpz_add(l2.im, temp1.im, l2.im);

    // 5. U' を計算．全体に (Z1 Z2)^6 がかかっている．
    // 18M, 1S
    ExtendedNumber f5Z2, f6U21, rV21;
    mpz_mul(f5Z2.re, this->f.coeff[5].re, m.Z.re);
    mpz_mod(f5Z2.re, f5Z2.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[5].im, m.Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(f5Z2.re, f5Z2.re, temp);
    mpz_mul(f5Z2.im, this->f.coeff[5].re, m.Z.im);
    mpz_mod(f5Z2.im, f5Z2.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[5].im, m.Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(f5Z2.im, f5Z2.im, temp);



    mpz_mul(f6U21.re, this->f.coeff[6].re, m.U1.re);
    mpz_mod(f6U21.re, f6U21.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[6].im, m.U1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(f6U21.re, f6U21.re, temp);
    mpz_mul(f6U21.im, this->f.coeff[6].re, m.U1.im);
    mpz_mod(f6U21.im, f6U21.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[6].im, m.U1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(f6U21.im, f6U21.im, temp);



    mpz_mul(rV21.re, r.re, m.V1.re);
    mpz_mod(rV21.re, rV21.re, ExtendedNumber::CHARA);
    mpz_mul(temp, r.im, m.V1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(rV21.re, rV21.re, temp);
    mpz_mul(rV21.im, r.re, m.V1.im);
    mpz_mod(rV21.im, rV21.im, ExtendedNumber::CHARA);
    mpz_mul(temp, r.im, m.V1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(rV21.im, rV21.im, temp);




    ExtendedNumber t4, t3, t2;
    mpz_mul(t4.re, rs.re, m.Z.re);
    mpz_mod(t4.re, t4.re, ExtendedNumber::CHARA);
    mpz_mul(temp, rs.im, m.Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(t4.re, t4.re, temp);
    mpz_mul(t4.im, rs.re, m.Z.im);
    mpz_mod(t4.im, t4.im, ExtendedNumber::CHARA);
    mpz_mul(temp, rs.im, m.Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(t4.im, t4.im, temp);



    mpz_mul(tempd.re, t4.re, this->f.coeff[6].re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, t4.im, this->f.coeff[6].im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, t4.re, this->f.coeff[6].im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, t4.im, this->f.coeff[6].re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(t4.re, tempd.re);
    mpz_set(t4.im, tempd.im);

    mpz_mul(temp1.re, s1.re, l3.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s1.im, l3.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, s1.re, l3.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s1.im, l3.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_sub(t4.re, temp1.re, t4.re);
    mpz_sub(t4.im, temp1.im, t4.im);
    mpz_mul(tempd.re, t4.re, m.Z.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, t4.im, m.Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, t4.re, m.Z.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, t4.im, m.Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(t4.re, tempd.re);
    mpz_set(t4.im, tempd.im);


    mpz_mul(t3.re, l2.re, s1.re);
    mpz_mod(t3.re, t3.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, s1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(t3.re, t3.re, temp);
    mpz_mul(t3.im, l2.re, s1.im);
    mpz_mod(t3.im, t3.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, s1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(t3.im, t3.im, temp);



    mpz_mul(temp1.re, l3.re, s0.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l3.im, s0.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, l3.re, s0.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l3.im, s0.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_add(t3.re, t3.re, temp1.re);
    mpz_add(t3.im, t3.im, temp1.im);
    mpz_sub(temp1.re, f5Z2.re, f6U21.re);
    mpz_sub(temp1.im, f5Z2.im, f6U21.im);
    mpz_mul(tempd.re, temp1.re, rs.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, rs.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, temp1.re, rs.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, rs.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(temp1.re, tempd.re);
    mpz_set(temp1.im, tempd.im);

    mpz_sub(t3.re, t3.re, temp1.re);
    mpz_sub(t3.im, t3.im, temp1.im);
    mpz_mul(tempd.re, t3.re, m.Z.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, t3.im, m.Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, t3.re, m.Z.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, t3.im, m.Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(t3.re, tempd.re);
    mpz_set(t3.im, tempd.im);


    mpz_add(t2.re, l1.re, rV21.re);
    mpz_add(t2.im, l1.im, rV21.im);
    mpz_add(t2.re, t2.re, rV21.re);
    mpz_add(t2.im, t2.im, rV21.im);
    mpz_mul(tempd.re, t2.re, s1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, t2.im, s1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, t2.re, s1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, t2.im, s1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(t2.re, tempd.re);
    mpz_set(t2.im, tempd.im);

    mpz_mul(temp1.re, s0.re, l2.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s0.im, l2.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, s0.re, l2.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s0.im, l2.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_add(t2.re, t2.re, temp1.re);
    mpz_add(t2.im, t2.im, temp1.im);
    mpz_mul(tempd.re, t2.re, m.Z.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, t2.im, m.Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, t2.re, m.Z.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, t2.im, m.Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(t2.re, tempd.re);
    mpz_set(t2.im, tempd.im);

    mpz_mul(temp1.re, this->f.coeff[4].re, m.Z.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[4].im, m.Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, this->f.coeff[4].re, m.Z.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[4].im, m.Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_mul(temp2.re, this->f.coeff[6].re, m.U0.re);
    mpz_mod(temp2.re, temp2.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[6].im, m.U0.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp2.re, temp2.re, temp);
    mpz_mul(temp2.im, this->f.coeff[6].re, m.U0.im);
    mpz_mod(temp2.im, temp2.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[6].im, m.U0.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp2.im, temp2.im, temp);



    mpz_sub(temp1.re, temp1.re, temp2.re);
    mpz_sub(temp1.im, temp1.im, temp2.im);
    mpz_mul(tempd.re, temp1.re, m.Z.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, m.Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, temp1.re, m.Z.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, m.Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(temp1.re, tempd.re);
    mpz_set(temp1.im, tempd.im);

    mpz_sub(temp2.re, f5Z2.re, f6U21.re);
    mpz_sub(temp2.im, f5Z2.im, f6U21.im);
    mpz_mul(tempd.re, m.U1.re, temp2.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, m.U1.im, temp2.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, m.U1.re, temp2.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, m.U1.im, temp2.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(temp2.re, tempd.re);
    mpz_set(temp2.im, tempd.im);

    mpz_sub(temp1.re, temp1.re, temp2.re);
    mpz_sub(temp1.im, temp1.im, temp2.im);
    mpz_mul(tempd.re, temp1.re, rs.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, rs.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, temp1.re, rs.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, rs.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(temp1.re, tempd.re);
    mpz_set(temp1.im, tempd.im);

    mpz_sub(t2.re, t2.re, temp1.re);
    mpz_sub(t2.im, t2.im, temp1.im);
    
    // 8M, 2S
    ExtendedNumber t4U11, t3Z1, Ud2, Ud1, Ud0, ZdS;
    mpz_mul(t4U11.re, t4.re, this->U1.re);
    mpz_mod(t4U11.re, t4U11.re, ExtendedNumber::CHARA);
    mpz_mul(temp, t4.im, this->U1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(t4U11.re, t4U11.re, temp);
    mpz_mul(t4U11.im, t4.re, this->U1.im);
    mpz_mod(t4U11.im, t4U11.im, ExtendedNumber::CHARA);
    mpz_mul(temp, t4.im, this->U1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(t4U11.im, t4U11.im, temp);



    mpz_mul(t3Z1.re, t3.re, this->Z.re);
    mpz_mod(t3Z1.re, t3Z1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, t3.im, this->Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(t3Z1.re, t3Z1.re, temp);
    mpz_mul(t3Z1.im, t3.re, this->Z.im);
    mpz_mod(t3Z1.im, t3Z1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, t3.im, this->Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(t3Z1.im, t3Z1.im, temp);



    mpz_mul(Ud2.re, t4.re, this->Z.re);
    mpz_mod(Ud2.re, Ud2.re, ExtendedNumber::CHARA);
    mpz_mul(temp, t4.im, this->Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(Ud2.re, Ud2.re, temp);
    mpz_mul(Ud2.im, t4.re, this->Z.im);
    mpz_mod(Ud2.im, Ud2.im, ExtendedNumber::CHARA);
    mpz_mul(temp, t4.im, this->Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(Ud2.im, Ud2.im, temp);



    mpz_mul(tempd.re, Ud2.re, this->Z.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud2.im, this->Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, Ud2.re, this->Z.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud2.im, this->Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(Ud2.re, tempd.re);
    mpz_set(Ud2.im, tempd.im);

    mpz_sub(Ud1.re, t3Z1.re, t4U11.re);
    mpz_sub(Ud1.im, t3Z1.im, t4U11.im);
    mpz_mul(tempd.re, Ud1.re, this->Z.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud1.im, this->Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, Ud1.re, this->Z.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud1.im, this->Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(Ud1.re, tempd.re);
    mpz_set(Ud1.im, tempd.im);

    mpz_mul(Ud0.re, t2.re, this->Z.re);
    mpz_mod(Ud0.re, Ud0.re, ExtendedNumber::CHARA);
    mpz_mul(temp, t2.im, this->Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(Ud0.re, Ud0.re, temp);
    mpz_mul(Ud0.im, t2.re, this->Z.im);
    mpz_mod(Ud0.im, Ud0.im, ExtendedNumber::CHARA);
    mpz_mul(temp, t2.im, this->Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(Ud0.im, Ud0.im, temp);



    mpz_mul(temp1.re, t4.re, this->U0.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, t4.im, this->U0.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, t4.re, this->U0.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, t4.im, this->U0.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_sub(Ud0.re, Ud0.re, temp1.re);
    mpz_sub(Ud0.im, Ud0.im, temp1.im);
    mpz_mul(tempd.re, Ud0.re, this->Z.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud0.im, this->Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, Ud0.re, this->Z.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud0.im, this->Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(Ud0.re, tempd.re);
    mpz_set(Ud0.im, tempd.im);
    mpz_sub(temp1.re, t3Z1.re, t4U11.re);
    mpz_sub(temp1.im, t3Z1.im, t4U11.im);
    mpz_mul(tempd.re, temp1.re, this->U1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, this->U1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, temp1.re, this->U1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, this->U1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(temp1.re, tempd.re);
    mpz_set(temp1.im, tempd.im);

    mpz_sub(Ud0.re, Ud0.re, temp1.re);
    mpz_sub(Ud0.im, Ud0.im, temp1.im);
    mpz_mul(ZdS.re, Ud2.re, Ud2.re);
    mpz_mod(ZdS.re, ZdS.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud2.im, Ud2.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(ZdS.re, ZdS.re, temp);
    mpz_mul(ZdS.im, Ud2.re, Ud2.im);
    mpz_mod(ZdS.im, ZdS.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud2.im, Ud2.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(ZdS.im, ZdS.im, temp);




    // 6. V' を計算．
    // 10M, 1S
    ExtendedNumber Vd0, Vd1;
    mpz_mul(Vd0.re, Ud1.re, l3.re);
    mpz_mod(Vd0.re, Vd0.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud1.im, l3.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(Vd0.re, Vd0.re, temp);
    mpz_mul(Vd0.im, Ud1.re, l3.im);
    mpz_mod(Vd0.im, Vd0.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud1.im, l3.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(Vd0.im, Vd0.im, temp);



    mpz_mul(temp1.re, l2.re, Ud2.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, Ud2.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, l2.re, Ud2.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, Ud2.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_sub(Vd0.re, Vd0.re, temp1.re);
    mpz_sub(Vd0.im, Vd0.im, temp1.im);
    mpz_mul(tempd.re, Vd0.re, Ud0.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Vd0.im, Ud0.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, Vd0.re, Ud0.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Vd0.im, Ud0.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(Vd0.re, tempd.re);
    mpz_set(Vd0.im, tempd.im);

    mpz_mul(temp1.re, m.V0.re, r.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, m.V0.im, r.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, m.V0.re, r.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, m.V0.im, r.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_add(temp1.re, temp1.re, l0.re);
    mpz_add(temp1.im, temp1.im, l0.im);
    mpz_mul(tempd.re, temp1.re, ZdS.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, ZdS.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, temp1.re, ZdS.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, ZdS.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(temp1.re, tempd.re);
    mpz_set(temp1.im, tempd.im);

    mpz_add(Vd0.re, Vd0.re, temp1.re);
    mpz_add(Vd0.im, Vd0.im, temp1.im);
    mpz_neg(Vd0.re, Vd0.re);
    mpz_neg(Vd0.im, Vd0.im);

    mpz_mul(Vd1.re, Ud1.re, Ud1.re);
    mpz_mod(Vd1.re, Vd1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud1.im, Ud1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(Vd1.re, Vd1.re, temp);
    mpz_mul(Vd1.im, Ud1.re, Ud1.im);
    mpz_mod(Vd1.im, Vd1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud1.im, Ud1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(Vd1.im, Vd1.im, temp);



    mpz_mul(temp1.re, Ud0.re, Ud2.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud0.im, Ud2.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, Ud0.re, Ud2.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud0.im, Ud2.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_sub(Vd1.re, Vd1.re, temp1.re);
    mpz_sub(Vd1.im, Vd1.im, temp1.im);
    mpz_mul(tempd.re, l3.re, Vd1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l3.im, Vd1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, l3.re, Vd1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l3.im, Vd1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(Vd1.re, tempd.re);
    mpz_set(Vd1.im, tempd.im);

    mpz_mul(temp1.re, Ud2.re, Ud1.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud2.im, Ud1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, Ud2.re, Ud1.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud2.im, Ud1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_mul(tempd.re, temp1.re, l2.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, l2.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, temp1.re, l2.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, l2.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(temp1.re, tempd.re);
    mpz_set(temp1.im, tempd.im);

    mpz_sub(Vd1.re, temp1.re, Vd1.re);
    mpz_sub(Vd1.im, temp1.im, Vd1.im);
    mpz_add(temp1.re, l1.re, rV21.re);
    mpz_add(temp1.im, l1.im, rV21.im);
    mpz_mul(tempd.re, temp1.re, ZdS.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, ZdS.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, temp1.re, ZdS.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, ZdS.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(temp1.re, tempd.re);
    mpz_set(temp1.im, tempd.im);

    mpz_sub(Vd1.re, Vd1.re, temp1.re);
    mpz_sub(Vd1.im, Vd1.im, temp1.im);

    // 7. Z' を調整．
    // 5M
    ExtendedNumber M, Zd;
    mpz_mul(M.re, Ud2.re, m.Z.re);
    mpz_mod(M.re, M.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud2.im, m.Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(M.re, M.re, temp);
    mpz_mul(M.im, Ud2.re, m.Z.im);
    mpz_mod(M.im, M.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud2.im, m.Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(M.im, M.im, temp);



    mpz_mul(tempd.re, M.re, r.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, M.im, r.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, M.re, r.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, M.im, r.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(M.re, tempd.re);
    mpz_set(M.im, tempd.im);

    mpz_mul(tempd.re, Ud1.re, M.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud1.im, M.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, Ud1.re, M.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud1.im, M.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(Ud1.re, tempd.re);
    mpz_set(Ud1.im, tempd.im);

    mpz_mul(tempd.re, Ud0.re, M.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud0.im, M.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, Ud0.re, M.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud0.im, M.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(Ud0.re, tempd.re);
    mpz_set(Ud0.im, tempd.im);

    mpz_mul(Zd.re, Ud2.re, M.re);
    mpz_mod(Zd.re, Zd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud2.im, M.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(Zd.re, Zd.re, temp);
    mpz_mul(Zd.im, Ud2.re, M.im);
    mpz_mod(Zd.im, Zd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud2.im, M.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(Zd.im, Zd.im, temp);




    ExtendedProjectiveMumford ret(f, h, Ud1, Ud0, Vd1, Vd0, Zd);
    return ret;
}

ExtendedProjectiveMumford ExtendedProjectiveMumford::LangeDoubling() const{
    //std::cout << "Projective Lange Doubling." << std::endl;
    // 59M, 9S

    mpz_t temp;
    mpz_init(temp);
    ExtendedNumber temp1, temp2, tempd;

    // 1. precomputation.
    // 3M, 2S

    ExtendedNumber U0Z, f5Z, Z2, Z3, Z4;

    mpz_mul(U0Z.re, this->U0.re, this->Z.re);
    mpz_mod(U0Z.re, U0Z.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->U0.im, this->Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(U0Z.re, U0Z.re, temp);
    mpz_mul(U0Z.im, this->U0.re, this->Z.im);
    mpz_mod(U0Z.im, U0Z.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->U0.im, this->Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(U0Z.im, U0Z.im, temp);



    mpz_mul(f5Z.re, this->f.coeff[5].re, this->Z.re);
    mpz_mod(f5Z.re, f5Z.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[5].im, this->Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(f5Z.re, f5Z.re, temp);
    mpz_mul(f5Z.im, this->f.coeff[5].re, this->Z.im);
    mpz_mod(f5Z.im, f5Z.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[5].im, this->Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(f5Z.im, f5Z.im, temp);



    mpz_mul(Z2.re, this->Z.re, this->Z.re);
    mpz_mod(Z2.re, Z2.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->Z.im, this->Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(Z2.re, Z2.re, temp);
    mpz_mul(Z2.im, this->Z.re, this->Z.im);
    mpz_mod(Z2.im, Z2.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->Z.im, this->Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(Z2.im, Z2.im, temp);



    mpz_mul(Z3.re, Z2.re, this->Z.re);
    mpz_mod(Z3.re, Z3.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Z2.im, this->Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(Z3.re, Z3.re, temp);
    mpz_mul(Z3.im, Z2.re, this->Z.im);
    mpz_mod(Z3.im, Z3.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Z2.im, this->Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(Z3.im, Z3.im, temp);



    mpz_mul(Z4.re, Z2.re, Z2.re);
    mpz_mod(Z4.re, Z4.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Z2.im, Z2.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(Z4.re, Z4.re, temp);
    mpz_mul(Z4.im, Z2.re, Z2.im);
    mpz_mod(Z4.im, Z4.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Z2.im, Z2.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(Z4.im, Z4.im, temp);




    // 2. v~ と u の終結式を計算．
    // 4M, 2S
    ExtendedNumber V1t, V0t;
    mpz_add(V1t.re, this->V1.re, this->V1.re);
    mpz_add(V1t.im, this->V1.im, this->V1.im);
    mpz_add(V0t.re, this->V0.re, this->V0.re);
    mpz_add(V0t.im, this->V0.im, this->V0.im);

    ExtendedNumber W0, W1, U1s, W2, W3, V0tZ, r;
    mpz_mul(W0.re, this->V1.re, this->V1.re);
    mpz_mod(W0.re, W0.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->V1.im, this->V1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(W0.re, W0.re, temp);
    mpz_mul(W0.im, this->V1.re, this->V1.im);
    mpz_mod(W0.im, W0.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->V1.im, this->V1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(W0.im, W0.im, temp);



    mpz_mul(U1s.re, this->U1.re, this->U1.re);
    mpz_mod(U1s.re, U1s.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->U1.im, this->U1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(U1s.re, U1s.re, temp);
    mpz_mul(U1s.im, this->U1.re, this->U1.im);
    mpz_mod(U1s.im, U1s.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->U1.im, this->U1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(U1s.im, U1s.im, temp);



    mpz_add(W2.re, W0.re, W0.re);
    mpz_add(W2.im, W0.im, W0.im);
    mpz_add(W2.re, W2.re, W2.re);
    mpz_add(W2.im, W2.im, W2.im);
    mpz_mul(W3.re, this->U1.re, V1t.re);
    mpz_mod(W3.re, W3.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->U1.im, V1t.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(W3.re, W3.re, temp);
    mpz_mul(W3.im, this->U1.re, V1t.im);
    mpz_mod(W3.im, W3.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->U1.im, V1t.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(W3.im, W3.im, temp);



    mpz_mul(V0tZ.re, V0t.re, this->Z.re);
    mpz_mod(V0tZ.re, V0tZ.re, ExtendedNumber::CHARA);
    mpz_mul(temp, V0t.im, this->Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(V0tZ.re, V0tZ.re, temp);
    mpz_mul(V0tZ.im, V0t.re, this->Z.im);
    mpz_mod(V0tZ.im, V0tZ.im, ExtendedNumber::CHARA);
    mpz_mul(temp, V0t.im, this->Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(V0tZ.im, V0tZ.im, temp);



    mpz_mul(r.re, this->U0.re, W2.re);
    mpz_mod(r.re, r.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->U0.im, W2.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(r.re, r.re, temp);
    mpz_mul(r.im, this->U0.re, W2.im);
    mpz_mod(r.im, r.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->U0.im, W2.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(r.im, r.im, temp);



    mpz_sub(temp1.re, V0tZ.re, W3.re);
    mpz_sub(temp1.im, V0tZ.im, W3.im);
    mpz_mul(tempd.re, temp1.re, V0t.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, V0t.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, temp1.re, V0t.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, V0t.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(temp1.re, tempd.re);
    mpz_set(temp1.im, tempd.im);

    mpz_add(r.re, r.re, temp1.re);
    mpz_add(r.im, r.im, temp1.im);

    // 3. almost inverse を計算．
    // Z がかかっている．
    //ExtendedNumber inv1d = -V1t;
    // Z^2 がかかっている．
    ExtendedNumber inv0d;
    mpz_sub(inv0d.re, V0tZ.re, W3.re);
    mpz_sub(inv0d.im, V0tZ.im, W3.im);

    // 4. k を計算．
    // 15M
    //ExtendedNumber k4 = f6;
    ExtendedNumber k4U1, k4U0Z, k3, k3U0Z, f4Z2, k2, k1, k0;
    mpz_mul(k4U1.re, this->f.coeff[6].re, this->U1.re);
    mpz_mod(k4U1.re, k4U1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[6].im, this->U1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(k4U1.re, k4U1.re, temp);
    mpz_mul(k4U1.im, this->f.coeff[6].re, this->U1.im);
    mpz_mod(k4U1.im, k4U1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[6].im, this->U1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(k4U1.im, k4U1.im, temp);



    mpz_mul(k4U0Z.re, this->f.coeff[6].re, U0Z.re);
    mpz_mod(k4U0Z.re, k4U0Z.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[6].im, U0Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(k4U0Z.re, k4U0Z.re, temp);
    mpz_mul(k4U0Z.im, this->f.coeff[6].re, U0Z.im);
    mpz_mod(k4U0Z.im, k4U0Z.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[6].im, U0Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(k4U0Z.im, k4U0Z.im, temp);



    mpz_sub(k3.re, f5Z.re, k4U1.re);
    mpz_sub(k3.im, f5Z.im, k4U1.im);
    mpz_mul(k3U0Z.re, k3.re, U0Z.re);
    mpz_mod(k3U0Z.re, k3U0Z.re, ExtendedNumber::CHARA);
    mpz_mul(temp, k3.im, U0Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(k3U0Z.re, k3U0Z.re, temp);
    mpz_mul(k3U0Z.im, k3.re, U0Z.im);
    mpz_mod(k3U0Z.im, k3U0Z.im, ExtendedNumber::CHARA);
    mpz_mul(temp, k3.im, U0Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(k3U0Z.im, k3U0Z.im, temp);



    mpz_mul(f4Z2.re, this->f.coeff[4].re, Z2.re);
    mpz_mod(f4Z2.re, f4Z2.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[4].im, Z2.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(f4Z2.re, f4Z2.re, temp);
    mpz_mul(f4Z2.im, this->f.coeff[4].re, Z2.im);
    mpz_mod(f4Z2.im, f4Z2.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[4].im, Z2.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(f4Z2.im, f4Z2.im, temp);



    mpz_mul(k2.re, k3.re, this->U1.re);
    mpz_mod(k2.re, k2.re, ExtendedNumber::CHARA);
    mpz_mul(temp, k3.im, this->U1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(k2.re, k2.re, temp);
    mpz_mul(k2.im, k3.re, this->U1.im);
    mpz_mod(k2.im, k2.im, ExtendedNumber::CHARA);
    mpz_mul(temp, k3.im, this->U1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(k2.im, k2.im, temp);



    mpz_sub(k2.re, f4Z2.re, k2.re);
    mpz_sub(k2.im, f4Z2.im, k2.im);
    mpz_sub(k2.re, k2.re, k4U0Z.re);
    mpz_sub(k2.im, k2.im, k4U0Z.im);

    mpz_mul(k1.re, this->f.coeff[3].re, Z3.re);
    mpz_mod(k1.re, k1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[3].im, Z3.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(k1.re, k1.re, temp);
    mpz_mul(k1.im, this->f.coeff[3].re, Z3.im);
    mpz_mod(k1.im, k1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[3].im, Z3.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(k1.im, k1.im, temp);



    mpz_sub(k1.re, k1.re, k3U0Z.re);
    mpz_sub(k1.im, k1.im, k3U0Z.im);
    mpz_mul(temp1.re, k2.re, this->U1.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, k2.im, this->U1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, k2.re, this->U1.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, k2.im, this->U1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_sub(k1.re, k1.re, temp1.re);
    mpz_sub(k1.im, k1.im, temp1.im);

    mpz_mul(k0.re, this->f.coeff[2].re, Z4.re);
    mpz_mod(k0.re, k0.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[2].im, Z4.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(k0.re, k0.re, temp);
    mpz_mul(k0.im, this->f.coeff[2].re, Z4.im);
    mpz_mod(k0.im, k0.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[2].im, Z4.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(k0.im, k0.im, temp);



    mpz_mul(temp1.re, W0.re, Z2.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, W0.im, Z2.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, W0.re, Z2.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, W0.im, Z2.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_sub(k0.re, k0.re, temp1.re);
    mpz_sub(k0.im, k0.im, temp1.im);
    mpz_mul(temp1.re, k2.re, U0Z.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, k2.im, U0Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, k2.re, U0Z.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, k2.im, U0Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_sub(k0.re, k0.re, temp1.re);
    mpz_sub(k0.im, k0.im, temp1.im);
    mpz_mul(temp1.re, k1.re, this->U1.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, k1.im, this->U1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, k1.re, this->U1.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, k1.im, this->U1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_sub(k0.re, k0.re, temp1.re);
    mpz_sub(k0.im, k0.im, temp1.im);

    ExtendedNumber k1d, k0d;
    mpz_add(k1d.re, k4U0Z.re, k4U0Z.re);
    mpz_add(k1d.im, k4U0Z.im, k4U0Z.im);
    mpz_sub(k1d.re, k1d.re, k2.re);
    mpz_sub(k1d.im, k1d.im, k2.im);
    mpz_mul(tempd.re, k1d.re, this->U1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, k1d.im, this->U1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, k1d.re, this->U1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, k1d.im, this->U1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(k1d.re, tempd.re);
    mpz_set(k1d.im, tempd.im);

    mpz_sub(temp1.re, k3.re, k4U1.re);
    mpz_sub(temp1.im, k3.im, k4U1.im);
    mpz_mul(tempd.re, temp1.re, U1s.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, U1s.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, temp1.re, U1s.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, U1s.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(temp1.re, tempd.re);
    mpz_set(temp1.im, tempd.im);

    mpz_add(k1d.re, k1d.re, temp1.re);
    mpz_add(k1d.im, k1d.im, temp1.im);
    mpz_add(k1d.re, k1.re, k1d.re);
    mpz_add(k1d.im, k1.im, k1d.im);
    mpz_sub(k1d.re, k1d.re, k3U0Z.re);
    mpz_sub(k1d.im, k1d.im, k3U0Z.im);

    // Z^4 がかかっている．
    mpz_sub(k0d.re, k3.re, k4U1.re);
    mpz_sub(k0d.im, k3.im, k4U1.im);
    mpz_mul(tempd.re, k0d.re, this->U1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, k0d.im, this->U1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, k0d.re, this->U1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, k0d.im, this->U1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(k0d.re, tempd.re);
    mpz_set(k0d.im, tempd.im);

    mpz_add(k0d.re, k0d.re, k4U0Z.re);
    mpz_add(k0d.im, k0d.im, k4U0Z.im);
    mpz_sub(k0d.re, k0d.re, k2.re);
    mpz_sub(k0d.im, k0d.im, k2.im);
    mpz_mul(tempd.re, k0d.re, U0Z.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, k0d.im, U0Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, k0d.re, U0Z.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, k0d.im, U0Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(k0d.re, tempd.re);
    mpz_set(k0d.im, tempd.im);

    mpz_add(k0d.re, k0d.re, k0.re);
    mpz_add(k0d.im, k0d.im, k0.im);

    // 5. s を計算．
    // 5M
    mpz_mul(W0.re, k0d.re, inv0d.re);
    mpz_mod(W0.re, W0.re, ExtendedNumber::CHARA);
    mpz_mul(temp, k0d.im, inv0d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(W0.re, W0.re, temp);
    mpz_mul(W0.im, k0d.re, inv0d.im);
    mpz_mod(W0.im, W0.im, ExtendedNumber::CHARA);
    mpz_mul(temp, k0d.im, inv0d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(W0.im, W0.im, temp);



    mpz_mul(W1.re, k1d.re, V1t.re);
    mpz_mod(W1.re, W1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, k1d.im, V1t.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(W1.re, W1.re, temp);
    mpz_mul(W1.im, k1d.re, V1t.im);
    mpz_mod(W1.im, W1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, k1d.im, V1t.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(W1.im, W1.im, temp);



    mpz_neg(W1.re, W1.re);
    mpz_neg(W1.im, W1.im);
    ExtendedNumber s1d, s0d;
    mpz_sub(s1d.re, inv0d.re, V1t.re);
    mpz_sub(s1d.im, inv0d.im, V1t.im);
    mpz_add(temp1.re, k0d.re, k1d.re);
    mpz_add(temp1.im, k0d.im, k1d.im);
    mpz_mul(tempd.re, s1d.re, temp1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, temp1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, s1d.re, temp1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, temp1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(s1d.re, tempd.re);
    mpz_set(s1d.im, tempd.im);

    mpz_sub(s1d.re, s1d.re, W0.re);
    mpz_sub(s1d.im, s1d.im, W0.im);
    mpz_add(temp1.re, ExtendedNumber::ONE().re, this->U1.re);
    mpz_add(temp1.im, ExtendedNumber::ONE().im, this->U1.im);
    mpz_mul(tempd.re, temp1.re, W1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, W1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, temp1.re, W1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, W1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(temp1.re, tempd.re);
    mpz_set(temp1.im, tempd.im);

    mpz_sub(s1d.re, s1d.re, temp1.re);
    mpz_sub(s1d.im, s1d.im, temp1.im);
    mpz_mul(s0d.re, W1.re, U0Z.re);
    mpz_mod(s0d.re, s0d.re, ExtendedNumber::CHARA);
    mpz_mul(temp, W1.im, U0Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(s0d.re, s0d.re, temp);
    mpz_mul(s0d.im, W1.re, U0Z.im);
    mpz_mod(s0d.im, s0d.im, ExtendedNumber::CHARA);
    mpz_mul(temp, W1.im, U0Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(s0d.im, s0d.im, temp);



    mpz_sub(s0d.re, W0.re, s0d.re);
    mpz_sub(s0d.im, W0.im, s0d.im);

    // 6. U' を計算．
    // 12M, 4S
    ExtendedNumber rs, f5Zk, Z3r, rsZ4;
    mpz_mul(rs.re, r.re, r.re);
    mpz_mod(rs.re, rs.re, ExtendedNumber::CHARA);
    mpz_mul(temp, r.im, r.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(rs.re, rs.re, temp);
    mpz_mul(rs.im, r.re, r.im);
    mpz_mod(rs.im, rs.im, ExtendedNumber::CHARA);
    mpz_mul(temp, r.im, r.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(rs.im, rs.im, temp);



    mpz_sub(f5Zk.re, f5Z.re, k4U1.re);
    mpz_sub(f5Zk.im, f5Z.im, k4U1.im);
    mpz_sub(f5Zk.re, f5Zk.re, k4U1.re);
    mpz_sub(f5Zk.im, f5Zk.im, k4U1.im);
    mpz_mul(Z3r.re, Z3.re, r.re);
    mpz_mod(Z3r.re, Z3r.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Z3.im, r.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(Z3r.re, Z3r.re, temp);
    mpz_mul(Z3r.im, Z3.re, r.im);
    mpz_mod(Z3r.im, Z3r.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Z3.im, r.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(Z3r.im, Z3r.im, temp);



    mpz_mul(rsZ4.re, rs.re, Z4.re);
    mpz_mod(rsZ4.re, rsZ4.re, ExtendedNumber::CHARA);
    mpz_mul(temp, rs.im, Z4.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(rsZ4.re, rsZ4.re, temp);
    mpz_mul(rsZ4.im, rs.re, Z4.im);
    mpz_mod(rsZ4.im, rsZ4.im, ExtendedNumber::CHARA);
    mpz_mul(temp, rs.im, Z4.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(rsZ4.im, rsZ4.im, temp);




    ExtendedNumber Ud2, Ud1, V1Z3r, Ud0, Zd, Zd2;
    mpz_mul(Ud2.re, s1d.re, s1d.re);
    mpz_mod(Ud2.re, Ud2.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, s1d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(Ud2.re, Ud2.re, temp);
    mpz_mul(Ud2.im, s1d.re, s1d.im);
    mpz_mod(Ud2.im, Ud2.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, s1d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(Ud2.im, Ud2.im, temp);



    mpz_mul(temp1.re, rsZ4.re, this->f.coeff[6].re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, rsZ4.im, this->f.coeff[6].im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, rsZ4.re, this->f.coeff[6].im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, rsZ4.im, this->f.coeff[6].re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_sub(Ud2.re, Ud2.re, temp1.re);
    mpz_sub(Ud2.im, Ud2.im, temp1.im);

    mpz_mul(Ud1.re, s1d.re, s0d.re);
    mpz_mod(Ud1.re, Ud1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, s0d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(Ud1.re, Ud1.re, temp);
    mpz_mul(Ud1.im, s1d.re, s0d.im);
    mpz_mod(Ud1.im, Ud1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, s0d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(Ud1.im, Ud1.im, temp);



    mpz_add(Ud1.re, Ud1.re, Ud1.re);
    mpz_add(Ud1.im, Ud1.im, Ud1.im);
    mpz_mul(temp1.re, rsZ4.re, f5Zk.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, rsZ4.im, f5Zk.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, rsZ4.re, f5Zk.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, rsZ4.im, f5Zk.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_sub(Ud1.re, Ud1.re, temp1.re);
    mpz_sub(Ud1.im, Ud1.im, temp1.im);

    mpz_mul(V1Z3r.re, this->V1.re, Z3r.re);
    mpz_mod(V1Z3r.re, V1Z3r.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->V1.im, Z3r.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(V1Z3r.re, V1Z3r.re, temp);
    mpz_mul(V1Z3r.im, this->V1.re, Z3r.im);
    mpz_mod(V1Z3r.im, V1Z3r.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->V1.im, Z3r.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(V1Z3r.im, V1Z3r.im, temp);




    mpz_add(Ud0.re, U0Z.re, U0Z.re);
    mpz_add(Ud0.im, U0Z.im, U0Z.im);
    mpz_add(Ud0.re, Ud0.re, U1s.re);
    mpz_add(Ud0.im, Ud0.im, U1s.im);
    mpz_mul(tempd.re, Ud0.re, this->f.coeff[6].re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud0.im, this->f.coeff[6].im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, Ud0.re, this->f.coeff[6].im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud0.im, this->f.coeff[6].re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(Ud0.re, tempd.re);
    mpz_set(Ud0.im, tempd.im);

    mpz_sub(Ud0.re, f4Z2.re, Ud0.re);
    mpz_sub(Ud0.im, f4Z2.im, Ud0.im);
    mpz_mul(temp1.re, this->U1.re, f5Zk.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->U1.im, f5Zk.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, this->U1.re, f5Zk.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->U1.im, f5Zk.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_add(temp1.re, temp1.re, temp1.re);
    mpz_add(temp1.im, temp1.im, temp1.im);
    mpz_sub(Ud0.re, Ud0.re, temp1.re);
    mpz_sub(Ud0.im, Ud0.im, temp1.im);
    mpz_mul(tempd.re, rsZ4.re, Ud0.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, rsZ4.im, Ud0.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, rsZ4.re, Ud0.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, rsZ4.im, Ud0.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(Ud0.re, tempd.re);
    mpz_set(Ud0.im, tempd.im);


    mpz_mul(temp1.re, s1d.re, V1Z3r.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, V1Z3r.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, s1d.re, V1Z3r.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, V1Z3r.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_add(temp1.re, temp1.re, temp1.re);
    mpz_add(temp1.im, temp1.im, temp1.im);
    mpz_sub(Ud0.re, temp1.re, Ud0.re);
    mpz_sub(Ud0.im, temp1.im, Ud0.im);

    mpz_mul(temp1.re, s0d.re, s0d.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s0d.im, s0d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, s0d.re, s0d.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s0d.im, s0d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_add(Ud0.re, Ud0.re, temp1.re);
    mpz_add(Ud0.im, Ud0.im, temp1.im);

    mpz_mul(tempd.re, Ud1.re, Z.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud1.im, Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, Ud1.re, Z.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud1.im, Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(Ud1.re, tempd.re);
    mpz_set(Ud1.im, tempd.im);

    mpz_mul(Zd.re, Ud2.re, Z2.re);
    mpz_mod(Zd.re, Zd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud2.im, Z2.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(Zd.re, Zd.re, temp);
    mpz_mul(Zd.im, Ud2.re, Z2.im);
    mpz_mod(Zd.im, Zd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud2.im, Z2.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(Zd.im, Zd.im, temp);



    mpz_mul(Zd2.re, Zd.re, Zd.re);
    mpz_mod(Zd2.re, Zd2.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Zd.im, Zd.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(Zd2.re, Zd2.re, temp);
    mpz_mul(Zd2.im, Zd.re, Zd.im);
    mpz_mod(Zd2.im, Zd2.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Zd.im, Zd.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(Zd2.im, Zd2.im, temp);




    // 7. V' を計算．
    // 15M, 1S
    // r Z^7 がかかっている．
    ExtendedNumber l3, l2, l1, l0, l2Zd;
    mpz_mul(l3.re, s1d.re, Z2.re);
    mpz_mod(l3.re, l3.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, Z2.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(l3.re, l3.re, temp);
    mpz_mul(l3.im, s1d.re, Z2.im);
    mpz_mod(l3.im, l3.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, Z2.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(l3.im, l3.im, temp);



    mpz_mul(l2.re, s1d.re, this->U1.re);
    mpz_mod(l2.re, l2.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, this->U1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(l2.re, l2.re, temp);
    mpz_mul(l2.im, s1d.re, this->U1.im);
    mpz_mod(l2.im, l2.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, this->U1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(l2.im, l2.im, temp);



    mpz_add(l2.re, l2.re, s0d.re);
    mpz_add(l2.im, l2.im, s0d.im);
    mpz_mul(tempd.re, l2.re, this->Z.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, this->Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, l2.re, this->Z.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, this->Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(l2.re, tempd.re);
    mpz_set(l2.im, tempd.im);

    mpz_mul(l1.re, s1d.re, U0Z.re);
    mpz_mod(l1.re, l1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, U0Z.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(l1.re, l1.re, temp);
    mpz_mul(l1.im, s1d.re, U0Z.im);
    mpz_mod(l1.im, l1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, U0Z.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(l1.im, l1.im, temp);



    mpz_mul(temp1.re, s0d.re, this->U1.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s0d.im, this->U1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, s0d.re, this->U1.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s0d.im, this->U1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_add(l1.re, l1.re, temp1.re);
    mpz_add(l1.im, l1.im, temp1.im);
    mpz_mul(l0.re, s0d.re, this->U0.re);
    mpz_mod(l0.re, l0.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s0d.im, this->U0.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(l0.re, l0.re, temp);
    mpz_mul(l0.im, s0d.re, this->U0.im);
    mpz_mod(l0.im, l0.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s0d.im, this->U0.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(l0.im, l0.im, temp);



    mpz_mul(l2Zd.re, l2.re, Zd.re);
    mpz_mod(l2Zd.re, l2Zd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, Zd.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(l2Zd.re, l2Zd.re, temp);
    mpz_mul(l2Zd.im, l2.re, Zd.im);
    mpz_mod(l2Zd.im, l2Zd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, Zd.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(l2Zd.im, l2Zd.im, temp);




    ExtendedNumber Vd1, Vd0;
    mpz_mul(Vd1.re, Ud1.re, Ud1.re);
    mpz_mod(Vd1.re, Vd1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud1.im, Ud1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(Vd1.re, Vd1.re, temp);
    mpz_mul(Vd1.im, Ud1.re, Ud1.im);
    mpz_mod(Vd1.im, Vd1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud1.im, Ud1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(Vd1.im, Vd1.im, temp);



    mpz_mul(temp1.re, Ud0.re, Zd.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud0.im, Zd.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, Ud0.re, Zd.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud0.im, Zd.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_sub(Vd1.re, Vd1.re, temp1.re);
    mpz_sub(Vd1.im, Vd1.im, temp1.im);
    mpz_mul(tempd.re, Vd1.re, l3.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Vd1.im, l3.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, Vd1.re, l3.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Vd1.im, l3.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(Vd1.re, tempd.re);
    mpz_set(Vd1.im, tempd.im);

    mpz_mul(temp1.re, Ud1.re, l2Zd.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud1.im, l2Zd.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, Ud1.re, l2Zd.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud1.im, l2Zd.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_sub(Vd1.re, temp1.re, Vd1.re);
    mpz_sub(Vd1.im, temp1.im, Vd1.im);
    mpz_add(temp1.re, l1.re, V1Z3r.re);
    mpz_add(temp1.im, l1.im, V1Z3r.im);
    mpz_mul(tempd.re, temp1.re, Zd2.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, Zd2.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, temp1.re, Zd2.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, Zd2.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(temp1.re, tempd.re);
    mpz_set(temp1.im, tempd.im);

    mpz_sub(Vd1.re, Vd1.re, temp1.re);
    mpz_sub(Vd1.im, Vd1.im, temp1.im);

    mpz_mul(Vd0.re, l3.re, Ud1.re);
    mpz_mod(Vd0.re, Vd0.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l3.im, Ud1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(Vd0.re, Vd0.re, temp);
    mpz_mul(Vd0.im, l3.re, Ud1.im);
    mpz_mod(Vd0.im, Vd0.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l3.im, Ud1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(Vd0.im, Vd0.im, temp);



    mpz_sub(Vd0.re, Vd0.re, l2Zd.re);
    mpz_sub(Vd0.im, Vd0.im, l2Zd.im);
    mpz_mul(tempd.re, Vd0.re, Ud0.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Vd0.im, Ud0.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, Vd0.re, Ud0.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Vd0.im, Ud0.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(Vd0.re, tempd.re);
    mpz_set(Vd0.im, tempd.im);

    mpz_mul(temp1.re, this->V0.re, Z3r.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->V0.im, Z3r.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, this->V0.re, Z3r.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->V0.im, Z3r.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_add(temp1.re, temp1.re, l0.re);
    mpz_add(temp1.im, temp1.im, l0.im);
    mpz_mul(tempd.re, temp1.re, Zd2.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, Zd2.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, temp1.re, Zd2.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, Zd2.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(temp1.re, tempd.re);
    mpz_set(temp1.im, tempd.im);

    mpz_add(Vd0.re, Vd0.re, temp1.re);
    mpz_add(Vd0.im, Vd0.im, temp1.im);
    mpz_neg(Vd0.re, Vd0.re);
    mpz_neg(Vd0.im, Vd0.im);

    // 8. 調整
    // 5M
    ExtendedNumber ZM;
    mpz_mul(ZM.re, Z4.re, Zd.re);
    mpz_mod(ZM.re, ZM.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Z4.im, Zd.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(ZM.re, ZM.re, temp);
    mpz_mul(ZM.im, Z4.re, Zd.im);
    mpz_mod(ZM.im, ZM.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Z4.im, Zd.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(ZM.im, ZM.im, temp);



    mpz_mul(tempd.re, ZM.re, r.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, ZM.im, r.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, ZM.re, r.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, ZM.im, r.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(ZM.re, tempd.re);
    mpz_set(ZM.im, tempd.im);

    mpz_mul(tempd.re, Ud1.re, ZM.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud1.im, ZM.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, Ud1.re, ZM.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud1.im, ZM.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(Ud1.re, tempd.re);
    mpz_set(Ud1.im, tempd.im);

    mpz_mul(tempd.re, Ud0.re, ZM.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud0.im, ZM.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, Ud0.re, ZM.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Ud0.im, ZM.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(Ud0.re, tempd.re);
    mpz_set(Ud0.im, tempd.im);

    mpz_mul(tempd.re, Zd.re, ZM.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, Zd.im, ZM.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, Zd.re, ZM.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, Zd.im, ZM.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(Zd.re, tempd.re);
    mpz_set(Zd.im, tempd.im);


    ExtendedProjectiveMumford ret(f, h, Ud1, Ud0, Vd1, Vd0, Zd);
    return ret;
}

ExtendedProjectiveMumford ExtendedProjectiveMumford::zero(const ExtendedPolynomial& f, const ExtendedPolynomial& h){
    ExtendedProjectiveMumford zero(f, h);
    return zero;
}

bool ExtendedProjectiveMumford::isZero() const{
    if(V1.isZero() && V0.isZero()){
        if(U1.isZero() && U0 == ExtendedNumber::ONE()){
            return true;
        }
    }
    return false;
}

void ExtendedProjectiveMumford::print() const{
    std::cerr << "[" << U1 << ", " << U0 << "]" << std::endl;
    std::cerr << "[" << V1 << ", " << V0 << "]" << std::endl;
    return;
}




