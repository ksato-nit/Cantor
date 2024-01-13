#include "extended_mumford.hpp"


ExtendedMumford::ExtendedMumford(ExtendedPolynomial f, ExtendedPolynomial h){
    this->f = f;
    this->h = h;
    this->u2 = ExtendedNumber::ONE();
    this->u1 = ExtendedNumber::ZERO();
    this->u0 = ExtendedNumber::ZERO();
    this->v1 = ExtendedNumber::ZERO();
    this->v0 = ExtendedNumber::ZERO();
    this->U1 = ExtendedNumber::ONE();
    this->U0 = ExtendedNumber::ZERO();
}

ExtendedMumford::ExtendedMumford(ExtendedPolynomial f, ExtendedPolynomial h, ExtendedNumber u1, ExtendedNumber u0, ExtendedNumber v1, ExtendedNumber v0){
    this->f = f;
    this->h = h;
    this->u2 = ExtendedNumber::ONE();
    this->u1 = u1;
    this->u0 = u0;
    this->v1 = v1;
    this->v0 = v0;
    this->U1 = u1 * u1;
    this->U0 = u1 * u0;
}

ExtendedMumford::ExtendedMumford(ExtendedPolynomial f, ExtendedPolynomial h, ExtendedNumber u1, ExtendedNumber u0, ExtendedNumber v1, ExtendedNumber v0, ExtendedNumber U1, ExtendedNumber U0){
    this->f = f;
    this->h = h;
    this->u2 = ExtendedNumber::ONE();
    this->u1 = u1;
    this->u0 = u0;
    this->v1 = v1;
    this->v0 = v0;
    this->U1 = U1;
    this->U0 = U0;
}

ExtendedMumford ExtendedMumford::CostelloScalarMultiple (const mpz_class& k_) const{
    ExtendedPolynomial f = this->f;
    ExtendedPolynomial h = this->h;
    mpz_class k = k_;

    // operator * の，Costello による計算．
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

    ExtendedMumford D = ExtendedMumford::zero(f, h);
    ExtendedMumford now = *this;
    for(int i = 0; i < count; ++i){
        if(bits[i] == 1){
            if(D.isZero()){
                D = now;
            }else if(now.isZero()){
                // D = D;
            }else{
                D = D.CostelloAdd(now);
            }
            //D.print();
        }
        now = now.CostelloDoubling();
    }
    return D;
}

ExtendedMumford ExtendedMumford::operator * (const mpz_class& k_) const{
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

    ExtendedMumford D = ExtendedMumford::zero(f, h);
    ExtendedMumford now = *this;
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

// todo: Mumford にあるのはおかしい
void ExtendedMumford::constant_invert(mpz_t op1, mpz_t op2, mpz_t p){
    // Fermat 法で逆元を計算する．
    // op1 = op2^(p-2) mod p;
    mpz_set_ui(op1, 1);

    mpz_t one;
    mpz_init_set_si(one, 1);

    mpz_t pow;
    mpz_init_set_si(pow, -2);
    mpz_add(pow, p, pow);

    mpz_t a;
    mpz_init_set(a, op2);

    mpz_t eval;
    mpz_init(eval);

    int n = mpz_sizeinbase(pow, 2);
    for(int i = 0; i <= n; ++i){
        mpz_and(eval, pow, one);
        if(mpz_sgn(eval) > 0){
            mpz_mul(op1, op1, a);
            mpz_mod(op1, op1, p);
        }
        mpz_mul(a, a, a);
        mpz_mod(a, a, p);
        mpz_fdiv_q_2exp(pow, pow, 1);        
    }
    return;
}


ExtendedMumford ExtendedMumford::CostelloAdd(const ExtendedMumford& m) const{
    //std::cerr << "Costello Addition." << std::endl;

    mpz_t temp;
    mpz_init(temp);
    ExtendedNumber temp1, tempd;

    ExtendedNumber u1S, v0D, v1D;
    mpz_add(u1S.re, this->u1.re, m.u1.re);
    mpz_add(u1S.im, this->u1.im, m.u1.im);
    mpz_sub(v0D.re, this->v0.re, m.v0.re);
    mpz_sub(v0D.im, this->v0.im, m.v0.im);
    mpz_sub(v1D.re, this->v1.re, m.v1.re);
    mpz_sub(v1D.im, this->v1.im, m.v1.im);

    ExtendedNumber M1, M2, M3, M4;
    mpz_sub(M1.re, this->U1.re, m.U1.re);
    mpz_sub(M1.im, this->U1.im, m.U1.im);
    mpz_sub(M1.re, M1.re, this->u0.re);
    mpz_sub(M1.im, M1.im, this->u0.im);
    mpz_add(M1.re, M1.re, m.u0.re);
    mpz_add(M1.im, M1.im, m.u0.im);
    mpz_sub(M2.re, m.U0.re, this->U0.re);
    mpz_sub(M2.im, m.U0.im, this->U0.im);
    mpz_sub(M3.re, this->u1.re, m.u1.re);
    mpz_sub(M3.im, this->u1.im, m.u1.im);
    mpz_sub(M4.re, m.u0.re, this->u0.re);
    mpz_sub(M4.im, m.u0.im, this->u0.im);

    ExtendedNumber t1, t2, t3, t4;
    mpz_sub(t1.re, M2.re, v0D.re);
    mpz_sub(t1.im, M2.im, v0D.im);
    mpz_sub(temp1.re, v1D.re, M1.re);
    mpz_sub(temp1.im, v1D.im, M1.im);
    mpz_mul(tempd.re, t1.re, temp1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, t1.im, temp1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, t1.re, temp1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, t1.im, temp1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(t1.re, tempd.re);
    mpz_set(t1.im, tempd.im);


    mpz_add(t2.re, M2.re, v0D.re);
    mpz_add(t2.im, M2.im, v0D.im);
    mpz_add(temp1.re, v1D.re, M1.re);
    mpz_add(temp1.im, v1D.im, M1.im);
    mpz_mul(tempd.re, t2.re, temp1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, t2.im, temp1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, t2.re, temp1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, t2.im, temp1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(t2.re, tempd.re);
    mpz_set(t2.im, tempd.im);

    mpz_neg(t2.re, t2.re);
    mpz_neg(t2.im, t2.im);

    mpz_sub(t3.re, M4.re, v0D.re);
    mpz_sub(t3.im, M4.im, v0D.im);
    mpz_sub(temp1.re, v1D.re, M3.re);
    mpz_sub(temp1.im, v1D.im, M3.im);
    mpz_mul(tempd.re, t3.re, temp1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, t3.im, temp1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, t3.re, temp1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, t3.im, temp1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(t3.re, tempd.re);
    mpz_set(t3.im, tempd.im);


    mpz_add(t4.re, M4.re, v0D.re);
    mpz_add(t4.im, M4.im, v0D.im);
    mpz_add(temp1.re, v1D.re, M3.re);
    mpz_add(temp1.im, v1D.im, M3.im);
    mpz_mul(tempd.re, t4.re, temp1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, t4.im, temp1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, t4.re, temp1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, t4.im, temp1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(t4.re, tempd.re);
    mpz_set(t4.im, tempd.im);

    mpz_neg(t4.re, t4.re);
    mpz_neg(t4.im, t4.im);

    ExtendedNumber l2_num, l3_num;
    mpz_sub(l2_num.re, t1.re, t2.re);
    mpz_sub(l2_num.im, t1.im, t2.im);
    mpz_sub(l3_num.re, t3.re, t4.re);
    mpz_sub(l3_num.im, t3.im, t4.im);

    ExtendedNumber d;
    mpz_sub(d.re, M4.re, M2.re);
    mpz_sub(d.im, M4.im, M2.im);
    mpz_add(temp1.re, M1.re, M3.re);
    mpz_add(temp1.im, M1.im, M3.im);
    mpz_mul(tempd.re, d.re, temp1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, d.im, temp1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, d.re, temp1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, d.im, temp1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(d.re, tempd.re);
    mpz_set(d.im, tempd.im);

    mpz_add(d.re, d.re, d.re);
    mpz_add(d.im, d.im, d.im);
    mpz_add(d.re, d.re, t3.re);
    mpz_add(d.im, d.im, t3.im);
    mpz_add(d.re, d.re, t4.re);
    mpz_add(d.im, d.im, t4.im);
    mpz_sub(d.re, d.re, t1.re);
    mpz_sub(d.im, d.im, t1.im);
    mpz_sub(d.re, d.re, t2.re);
    mpz_sub(d.im, d.im, t2.im);

    ExtendedNumber A, B, C, d_inv, d_shifted_inv;
    mpz_mul(A.re, d.re, d.re);
    mpz_mod(A.re, A.re, ExtendedNumber::CHARA);
    mpz_mul(temp, d.im, d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(A.re, A.re, temp);
    mpz_mul(A.im, d.re, d.im);
    mpz_mod(A.im, A.im, ExtendedNumber::CHARA);
    mpz_mul(temp, d.im, d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(A.im, A.im, temp);



    mpz_mul(temp1.re, this->f.coeff[6].re, A.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[6].im, A.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, this->f.coeff[6].re, A.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[6].im, A.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_mul(B.re, l3_num.re, l3_num.re);
    mpz_mod(B.re, B.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l3_num.im, l3_num.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(B.re, B.re, temp);
    mpz_mul(B.im, l3_num.re, l3_num.im);
    mpz_mod(B.im, B.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l3_num.im, l3_num.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(B.im, B.im, temp);



    mpz_sub(B.re, B.re, temp1.re);
    mpz_sub(B.im, B.im, temp1.im);
    mpz_mul(C.re, d.re, B.re);
    mpz_mod(C.re, C.re, ExtendedNumber::CHARA);
    mpz_mul(temp, d.im, B.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(C.re, C.re, temp);
    mpz_mul(C.im, d.re, B.im);
    mpz_mod(C.im, C.im, ExtendedNumber::CHARA);
    mpz_mul(temp, d.im, B.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(C.im, C.im, temp);



    mpz_t denom;
    mpz_init(denom);
    mpz_mul(temp, C.re, C.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_mul(denom, C.im, C.im);
    mpz_mod(denom, denom, ExtendedNumber::CHARA);
    mpz_add(denom, denom, temp);
    constant_invert(denom, denom, ExtendedNumber::CHARA);
    mpz_mul(C.re, C.re, denom);
    mpz_mod(C.re, C.re, ExtendedNumber::CHARA);
    mpz_mul(C.im, C.im, denom);
    mpz_mod(C.im, C.im, ExtendedNumber::CHARA);
    mpz_neg(C.im, C.im);

    mpz_mul(d_inv.re, B.re, C.re);
    mpz_mod(d_inv.re, d_inv.re, ExtendedNumber::CHARA);
    mpz_mul(temp, B.im, C.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(d_inv.re, d_inv.re, temp);
    mpz_mul(d_inv.im, B.re, C.im);
    mpz_mod(d_inv.im, d_inv.im, ExtendedNumber::CHARA);
    mpz_mul(temp, B.im, C.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(d_inv.im, d_inv.im, temp);



    mpz_mul(d_shifted_inv.re, d.re, A.re);
    mpz_mod(d_shifted_inv.re, d_shifted_inv.re, ExtendedNumber::CHARA);
    mpz_mul(temp, d.im, A.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(d_shifted_inv.re, d_shifted_inv.re, temp);
    mpz_mul(d_shifted_inv.im, d.re, A.im);
    mpz_mod(d_shifted_inv.im, d_shifted_inv.im, ExtendedNumber::CHARA);
    mpz_mul(temp, d.im, A.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(d_shifted_inv.im, d_shifted_inv.im, temp);



    mpz_mul(tempd.re, d_shifted_inv.re, C.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, d_shifted_inv.im, C.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, d_shifted_inv.re, C.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, d_shifted_inv.im, C.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(d_shifted_inv.re, tempd.re);
    mpz_set(d_shifted_inv.im, tempd.im);


    ExtendedNumber l2, l3;
    mpz_mul(l2.re, l2_num.re, d_inv.re);
    mpz_mod(l2.re, l2.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l2_num.im, d_inv.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(l2.re, l2.re, temp);
    mpz_mul(l2.im, l2_num.re, d_inv.im);
    mpz_mod(l2.im, l2.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l2_num.im, d_inv.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(l2.im, l2.im, temp);



    mpz_mul(l3.re, l3_num.re, d_inv.re);
    mpz_mod(l3.re, l3.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l3_num.im, d_inv.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(l3.re, l3.re, temp);
    mpz_mul(l3.im, l3_num.re, d_inv.im);
    mpz_mod(l3.im, l3.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l3_num.im, d_inv.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(l3.im, l3.im, temp);




    mpz_mul(temp1.re, l2.re, l3.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, l3.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, l2.re, l3.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, l3.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);




    ExtendedNumber u1dd;
    mpz_sub(u1dd.re, this->f.coeff[5].re, temp1.re);
    mpz_sub(u1dd.im, this->f.coeff[5].im, temp1.im);
    mpz_sub(u1dd.re, u1dd.re, temp1.re);
    mpz_sub(u1dd.im, u1dd.im, temp1.im);
    mpz_mul(tempd.re, u1dd.re, d_shifted_inv.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u1dd.im, d_shifted_inv.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, u1dd.re, d_shifted_inv.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u1dd.im, d_shifted_inv.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(u1dd.re, tempd.re);
    mpz_set(u1dd.im, tempd.im);

    mpz_add(u1dd.re, u1dd.re, u1S.re);
    mpz_add(u1dd.im, u1dd.im, u1S.im);
    mpz_neg(u1dd.re, u1dd.re);
    mpz_neg(u1dd.im, u1dd.im);

    ExtendedNumber u0dd;
    mpz_sub(u0dd.re, this->u0.re, this->U1.re);
    mpz_sub(u0dd.im, this->u0.im, this->U1.im);
    mpz_mul(tempd.re, u0dd.re, l3.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u0dd.im, l3.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, u0dd.re, l3.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u0dd.im, l3.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(u0dd.re, tempd.re);
    mpz_set(u0dd.im, tempd.im);

    mpz_mul(temp1.re, l2.re, this->u1.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, this->u1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, l2.re, this->u1.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, this->u1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_add(u0dd.re, u0dd.re, temp1.re);
    mpz_add(u0dd.im, u0dd.im, temp1.im);
    mpz_add(u0dd.re, u0dd.re, this->v1.re);
    mpz_add(u0dd.im, u0dd.im, this->v1.im);
    mpz_mul(tempd.re, u0dd.re, l3.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u0dd.im, l3.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, u0dd.re, l3.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u0dd.im, l3.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(u0dd.re, tempd.re);
    mpz_set(u0dd.im, tempd.im);

    mpz_add(u0dd.re, u0dd.re, u0dd.re);
    mpz_add(u0dd.im, u0dd.im, u0dd.im);
    mpz_mul(temp1.re, l2.re, l2.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, l2.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, l2.re, l2.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, l2.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_add(u0dd.re, u0dd.re, temp1.re);
    mpz_add(u0dd.im, u0dd.im, temp1.im);
    mpz_sub(u0dd.re, u0dd.re, this->f.coeff[4].re);
    mpz_sub(u0dd.im, u0dd.im, this->f.coeff[4].im);
    mpz_mul(tempd.re, u0dd.re, d_shifted_inv.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u0dd.im, d_shifted_inv.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, u0dd.re, d_shifted_inv.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u0dd.im, d_shifted_inv.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(u0dd.re, tempd.re);
    mpz_set(u0dd.im, tempd.im);

    mpz_mul(temp1.re, m.u1.re, this->u1.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, m.u1.im, this->u1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, m.u1.re, this->u1.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, m.u1.im, this->u1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_sub(u0dd.re, u0dd.re, temp1.re);
    mpz_sub(u0dd.im, u0dd.im, temp1.im);
    mpz_sub(u0dd.re, u0dd.re, this->u0.re);
    mpz_sub(u0dd.im, u0dd.im, this->u0.im);
    mpz_sub(u0dd.re, u0dd.re, m.u0.re);
    mpz_sub(u0dd.im, u0dd.im, m.u0.im);
    mpz_mul(temp1.re, u1S.re, u1dd.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u1S.im, u1dd.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, u1S.re, u1dd.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u1S.im, u1dd.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_sub(u0dd.re, u0dd.re, temp1.re);
    mpz_sub(u0dd.im, u0dd.im, temp1.im);

    ExtendedNumber U1dd, U0dd;
    mpz_mul(U1dd.re, u1dd.re, u1dd.re);
    mpz_mod(U1dd.re, U1dd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u1dd.im, u1dd.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(U1dd.re, U1dd.re, temp);
    mpz_mul(U1dd.im, u1dd.re, u1dd.im);
    mpz_mod(U1dd.im, U1dd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u1dd.im, u1dd.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(U1dd.im, U1dd.im, temp);



    mpz_mul(U0dd.re, u1dd.re, u0dd.re);
    mpz_mod(U0dd.re, U0dd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u1dd.im, u0dd.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(U0dd.re, U0dd.re, temp);
    mpz_mul(U0dd.im, u1dd.re, u0dd.im);
    mpz_mod(U0dd.im, U0dd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u1dd.im, u0dd.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(U0dd.im, U0dd.im, temp);




    ExtendedNumber v1dd, v0dd;
    mpz_sub(temp1.re, u1dd.re, this->u1.re);
    mpz_sub(temp1.im, u1dd.im, this->u1.im);
    mpz_mul(tempd.re, l2.re, temp1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, temp1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, l2.re, temp1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, temp1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(temp1.re, tempd.re);
    mpz_set(temp1.im, tempd.im);


    mpz_sub(v1dd.re, u0dd.re, U1dd.re);
    mpz_sub(v1dd.im, u0dd.im, U1dd.im);
    mpz_add(v1dd.re, v1dd.re, this->U1.re);
    mpz_add(v1dd.im, v1dd.im, this->U1.im);
    mpz_sub(v1dd.re, v1dd.re, this->u0.re);
    mpz_sub(v1dd.im, v1dd.im, this->u0.im);
    mpz_mul(tempd.re, v1dd.re, l3.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, v1dd.im, l3.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, v1dd.re, l3.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, v1dd.im, l3.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(v1dd.re, tempd.re);
    mpz_set(v1dd.im, tempd.im);

    mpz_add(v1dd.re, v1dd.re, temp1.re);
    mpz_add(v1dd.im, v1dd.im, temp1.im);
    mpz_sub(v1dd.re, v1dd.re, this->v1.re);
    mpz_sub(v1dd.im, v1dd.im, this->v1.im);

    mpz_sub(temp1.re, u0dd.re, this->u0.re);
    mpz_sub(temp1.im, u0dd.im, this->u0.im);
    mpz_mul(tempd.re, l2.re, temp1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, temp1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, l2.re, temp1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, temp1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(temp1.re, tempd.re);
    mpz_set(temp1.im, tempd.im);


    mpz_sub(v0dd.re, this->U0.re, U0dd.re);
    mpz_sub(v0dd.im, this->U0.im, U0dd.im);
    mpz_mul(tempd.re, v0dd.re, l3.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, v0dd.im, l3.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, v0dd.re, l3.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, v0dd.im, l3.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(v0dd.re, tempd.re);
    mpz_set(v0dd.im, tempd.im);

    mpz_add(v0dd.re, v0dd.re, temp1.re);
    mpz_add(v0dd.im, v0dd.im, temp1.im);
    mpz_sub(v0dd.re, v0dd.re, this->v0.re);
    mpz_sub(v0dd.im, v0dd.im, this->v0.im);

    return ExtendedMumford(f, h, u1dd, u0dd, v1dd, v0dd, U1dd, U0dd);
}


ExtendedMumford ExtendedMumford::LangeAdd(const ExtendedMumford& m) const{
    //std::cerr << "Lange Addition." << std::endl;

    mpz_t temp;
    mpz_init(temp);
    ExtendedNumber temp1, temp2, temp3, tempd;

    // 1. u1, u2 の終結式を計算．
    ExtendedNumber z1, z2, z3, r;
    mpz_sub(z1.re, this->u1.re, m.u1.re);
    mpz_sub(z1.im, this->u1.im, m.u1.im);
    mpz_sub(z2.re, m.u0.re, this->u0.re);
    mpz_sub(z2.im, m.u0.im, this->u0.im);
    mpz_mul(z3.re, this->u1.re, z1.re);
    mpz_mod(z3.re, z3.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->u1.im, z1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(z3.re, z3.re, temp);
    mpz_mul(z3.im, this->u1.re, z1.im);
    mpz_mod(z3.im, z3.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->u1.im, z1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(z3.im, z3.im, temp);



    mpz_add(z3.re, z3.re, z2.re);
    mpz_add(z3.im, z3.im, z2.im);
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



    mpz_mul(tempd.re, r.re, this->u0.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, r.im, this->u0.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, r.re, this->u0.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, r.im, this->u0.re);
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



    mpz_add(r.re, r.re, temp1.re);
    mpz_add(r.im, r.im, temp1.im);

    // 2. u2 の almost inverse (mod u1) を計算．
    //ExtendedNumber inv1 = z1;
    //ExtendedNumber inv0 = z3;

    // 3. s' を計算．
    ExtendedNumber w0, w1, w2, w3, s1d, s0d;
    mpz_sub(w0.re, this->v0.re, m.v0.re);
    mpz_sub(w0.im, this->v0.im, m.v0.im);
    mpz_sub(w1.re, this->v1.re, m.v1.re);
    mpz_sub(w1.im, this->v1.im, m.v1.im);
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




    mpz_add(temp1.re, w0.re, w1.re);
    mpz_add(temp1.im, w0.im, w1.im);
    mpz_add(s1d.re, z3.re, z1.re);
    mpz_add(s1d.im, z3.im, z1.im);
    mpz_mul(tempd.re, temp1.re, s1d.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, s1d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, temp1.re, s1d.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, s1d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(s1d.re, tempd.re);
    mpz_set(s1d.im, tempd.im);

    mpz_sub(s1d.re, s1d.re, w2.re);
    mpz_sub(s1d.im, s1d.im, w2.im);
    mpz_add(temp1.re, ExtendedNumber::ONE().re, this->u1.re);
    mpz_add(temp1.im, ExtendedNumber::ONE().im, this->u1.im);
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

    mpz_sub(s1d.re, s1d.re, temp1.re);
    mpz_sub(s1d.im, s1d.im, temp1.im);

    mpz_mul(s0d.re, this->u0.re, w3.re);
    mpz_mod(s0d.re, s0d.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->u0.im, w3.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(s0d.re, s0d.re, temp);
    mpz_mul(s0d.im, this->u0.re, w3.im);
    mpz_mod(s0d.im, s0d.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->u0.im, w3.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(s0d.im, s0d.im, temp);



    mpz_sub(s0d.re, w2.re, s0d.re);
    mpz_sub(s0d.im, w2.im, s0d.im);

    if(s1d.isZero()){
        //std::cerr << "Special case." << std::endl;
        // todo: ここの場合分けを厳密に書く．
        ExtendedMumford ret(f, h);
        return ret;
    }

    // 4. l' を計算．
    //ExtendedNumber l3d = s1d;
    ExtendedNumber l2d, l0d, l1d;
    mpz_mul(l2d.re, m.u1.re, s1d.re);
    mpz_mod(l2d.re, l2d.re, ExtendedNumber::CHARA);
    mpz_mul(temp, m.u1.im, s1d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(l2d.re, l2d.re, temp);
    mpz_mul(l2d.im, m.u1.re, s1d.im);
    mpz_mod(l2d.im, l2d.im, ExtendedNumber::CHARA);
    mpz_mul(temp, m.u1.im, s1d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(l2d.im, l2d.im, temp);



    mpz_mul(l0d.re, m.u0.re, s0d.re);
    mpz_mod(l0d.re, l0d.re, ExtendedNumber::CHARA);
    mpz_mul(temp, m.u0.im, s0d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(l0d.re, l0d.re, temp);
    mpz_mul(l0d.im, m.u0.re, s0d.im);
    mpz_mod(l0d.im, l0d.im, ExtendedNumber::CHARA);
    mpz_mul(temp, m.u0.im, s0d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(l0d.im, l0d.im, temp);



    mpz_add(l1d.re, s1d.re, s0d.re);
    mpz_add(l1d.im, s1d.im, s0d.im);
    mpz_add(temp1.re, m.u1.re, m.u0.re);
    mpz_add(temp1.im, m.u1.im, m.u0.im);
    mpz_mul(tempd.re, l1d.re, temp1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l1d.im, temp1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, l1d.re, temp1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l1d.im, temp1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(l1d.re, tempd.re);
    mpz_set(l1d.im, tempd.im);

    mpz_sub(l1d.re, l1d.re, l2d.re);
    mpz_sub(l1d.im, l1d.im, l2d.im);
    mpz_sub(l1d.re, l1d.re, l0d.re);
    mpz_sub(l1d.im, l1d.im, l0d.im);
    mpz_add(l2d.re, l2d.re, s0d.re);
    mpz_add(l2d.im, l2d.im, s0d.im);

    // 5. u' を計算．
    //ExtendedNumber k4 = f6;
    ExtendedNumber k3, k2;
    mpz_mul(k3.re, this->f.coeff[6].re, m.u1.re);
    mpz_mod(k3.re, k3.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[6].im, m.u1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(k3.re, k3.re, temp);
    mpz_mul(k3.im, this->f.coeff[6].re, m.u1.im);
    mpz_mod(k3.im, k3.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[6].im, m.u1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(k3.im, k3.im, temp);



    mpz_sub(k3.re, this->f.coeff[5].re, k3.re);
    mpz_sub(k3.im, this->f.coeff[5].im, k3.im);
    mpz_mul(k2.re, this->f.coeff[6].re, m.u0.re);
    mpz_mod(k2.re, k2.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[6].im, m.u0.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(k2.re, k2.re, temp);
    mpz_mul(k2.im, this->f.coeff[6].re, m.u0.im);
    mpz_mod(k2.im, k2.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[6].im, m.u0.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(k2.im, k2.im, temp);



    mpz_mul(temp1.re, k3.re, m.u1.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, k3.im, m.u1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, k3.re, m.u1.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, k3.im, m.u1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);


    mpz_sub(k2.re, this->f.coeff[4].re, k2.re);
    mpz_sub(k2.im, this->f.coeff[4].im, k2.im);
    mpz_sub(k2.re, k2.re, temp1.re);
    mpz_sub(k2.im, k2.im, temp1.im);

    ExtendedNumber t4, t3, t2;
    mpz_mul(temp1.re, r.re, r.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, r.im, r.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, r.re, r.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, r.im, r.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_mul(temp2.re, temp1.re, this->f.coeff[6].re);
    mpz_mod(temp2.re, temp2.re, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, this->f.coeff[6].im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp2.re, temp2.re, temp);
    mpz_mul(temp2.im, temp1.re, this->f.coeff[6].im);
    mpz_mod(temp2.im, temp2.im, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, this->f.coeff[6].re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp2.im, temp2.im, temp);



    mpz_mul(t4.re, s1d.re, s1d.re);
    mpz_mod(t4.re, t4.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, s1d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(t4.re, t4.re, temp);
    mpz_mul(t4.im, s1d.re, s1d.im);
    mpz_mod(t4.im, t4.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, s1d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(t4.im, t4.im, temp);



    mpz_sub(t4.re, t4.re, temp2.re);
    mpz_sub(t4.im, t4.im, temp2.im);
    mpz_mul(t3.re, s1d.re, l2d.re);
    mpz_mod(t3.re, t3.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, l2d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(t3.re, t3.re, temp);
    mpz_mul(t3.im, s1d.re, l2d.im);
    mpz_mod(t3.im, t3.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, l2d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(t3.im, t3.im, temp);



    mpz_mul(temp2.re, s0d.re, s1d.re);
    mpz_mod(temp2.re, temp2.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s0d.im, s1d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp2.re, temp2.re, temp);
    mpz_mul(temp2.im, s0d.re, s1d.im);
    mpz_mod(temp2.im, temp2.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s0d.im, s1d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp2.im, temp2.im, temp);



    mpz_add(t3.re, t3.re, temp2.re);
    mpz_add(t3.im, t3.im, temp2.im);
    mpz_mul(temp2.re, temp1.re, k3.re);
    mpz_mod(temp2.re, temp2.re, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, k3.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp2.re, temp2.re, temp);
    mpz_mul(temp2.im, temp1.re, k3.im);
    mpz_mod(temp2.im, temp2.im, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, k3.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp2.im, temp2.im, temp);



    mpz_sub(t3.re, t3.re, temp2.re);
    mpz_sub(t3.im, t3.im, temp2.im);
    mpz_mul(t2.re, r.re, m.v1.re);
    mpz_mod(t2.re, t2.re, ExtendedNumber::CHARA);
    mpz_mul(temp, r.im, m.v1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(t2.re, t2.re, temp);
    mpz_mul(t2.im, r.re, m.v1.im);
    mpz_mod(t2.im, t2.im, ExtendedNumber::CHARA);
    mpz_mul(temp, r.im, m.v1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(t2.im, t2.im, temp);



    mpz_add(t2.re, t2.re, t2.re);
    mpz_add(t2.im, t2.im, t2.im);
    mpz_add(t2.re, t2.re, l1d.re);
    mpz_add(t2.im, t2.im, l1d.im);
    mpz_mul(tempd.re, t2.re, s1d.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, t2.im, s1d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, t2.re, s1d.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, t2.im, s1d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(t2.re, tempd.re);
    mpz_set(t2.im, tempd.im);

    mpz_mul(temp2.re, s0d.re, l2d.re);
    mpz_mod(temp2.re, temp2.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s0d.im, l2d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp2.re, temp2.re, temp);
    mpz_mul(temp2.im, s0d.re, l2d.im);
    mpz_mod(temp2.im, temp2.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s0d.im, l2d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp2.im, temp2.im, temp);



    mpz_add(t2.re, t2.re, temp2.re);
    mpz_add(t2.im, t2.im, temp2.im);
    mpz_mul(temp2.re, temp1.re, k2.re);
    mpz_mod(temp2.re, temp2.re, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, k2.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp2.re, temp2.re, temp);
    mpz_mul(temp2.im, temp1.re, k2.im);
    mpz_mod(temp2.im, temp2.im, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, k2.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp2.im, temp2.im, temp);



    mpz_sub(t2.re, t2.re, temp2.re);
    mpz_sub(t2.im, t2.im, temp2.im);

    //ExtendedNumber u2d = t4;
    ExtendedNumber u1d, u0d;
    mpz_mul(u1d.re, t4.re, this->u1.re);
    mpz_mod(u1d.re, u1d.re, ExtendedNumber::CHARA);
    mpz_mul(temp, t4.im, this->u1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(u1d.re, u1d.re, temp);
    mpz_mul(u1d.im, t4.re, this->u1.im);
    mpz_mod(u1d.im, u1d.im, ExtendedNumber::CHARA);
    mpz_mul(temp, t4.im, this->u1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(u1d.im, u1d.im, temp);



    mpz_sub(u1d.re, t3.re, u1d.re);
    mpz_sub(u1d.im, t3.im, u1d.im);
    mpz_mul(u0d.re, u1d.re, this->u1.re);
    mpz_mod(u0d.re, u0d.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u1d.im, this->u1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(u0d.re, u0d.re, temp);
    mpz_mul(u0d.im, u1d.re, this->u1.im);
    mpz_mod(u0d.im, u0d.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u1d.im, this->u1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(u0d.im, u0d.im, temp);



    mpz_mul(temp1.re, t4.re, this->u0.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, t4.im, this->u0.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, t4.re, this->u0.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, t4.im, this->u0.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_add(u0d.re, u0d.re, temp1.re);
    mpz_add(u0d.im, u0d.im, temp1.im);
    mpz_sub(u0d.re, t2.re, u0d.re);
    mpz_sub(u0d.im, t2.im, u0d.im);

    // 6. u' をモニックにする．
    mpz_mul(w1.re, t4.re, r.re);
    mpz_mod(w1.re, w1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, t4.im, r.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(w1.re, w1.re, temp);
    mpz_mul(w1.im, t4.re, r.im);
    mpz_mod(w1.im, w1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, t4.im, r.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(w1.im, w1.im, temp);



    mpz_t denom;
    mpz_init(denom);
    mpz_mul(temp, w1.re, w1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_mul(denom, w1.im, w1.im);
    mpz_mod(denom, denom, ExtendedNumber::CHARA);
    mpz_add(denom, denom, temp);
    constant_invert(denom, denom, ExtendedNumber::CHARA);
    mpz_mul(w1.re, w1.re, denom);
    mpz_mod(w1.re, w1.re, ExtendedNumber::CHARA);
    mpz_mul(w1.im, w1.im, denom);
    mpz_mod(w1.im, w1.im, ExtendedNumber::CHARA);
    mpz_neg(w1.im, w1.im);

    mpz_mul(w2.re, w1.re, r.re);
    mpz_mod(w2.re, w2.re, ExtendedNumber::CHARA);
    mpz_mul(temp, w1.im, r.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(w2.re, w2.re, temp);
    mpz_mul(w2.im, w1.re, r.im);
    mpz_mod(w2.im, w2.im, ExtendedNumber::CHARA);
    mpz_mul(temp, w1.im, r.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(w2.im, w2.im, temp);



    mpz_mul(w3.re, w1.re, t4.re);
    mpz_mod(w3.re, w3.re, ExtendedNumber::CHARA);
    mpz_mul(temp, w1.im, t4.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(w3.re, w3.re, temp);
    mpz_mul(w3.im, w1.re, t4.im);
    mpz_mod(w3.im, w3.im, ExtendedNumber::CHARA);
    mpz_mul(temp, w1.im, t4.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(w3.im, w3.im, temp);



    mpz_mul(tempd.re, u1d.re, w2.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u1d.im, w2.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, u1d.re, w2.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u1d.im, w2.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(u1d.re, tempd.re);
    mpz_set(u1d.im, tempd.im);

    mpz_mul(tempd.re, u0d.re, w2.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u0d.im, w2.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, u0d.re, w2.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u0d.im, w2.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(u0d.re, tempd.re);
    mpz_set(u0d.im, tempd.im);


    // 7. v' を計算．
    ExtendedNumber v1d, v0d;

    mpz_mul(v1d.re, u1d.re, u1d.re);
    mpz_mod(v1d.re, v1d.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u1d.im, u1d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(v1d.re, v1d.re, temp);
    mpz_mul(v1d.im, u1d.re, u1d.im);
    mpz_mod(v1d.im, v1d.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u1d.im, u1d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(v1d.im, v1d.im, temp);



    mpz_sub(v1d.re, u0d.re, v1d.re);
    mpz_sub(v1d.im, u0d.im, v1d.im);
    mpz_mul(tempd.re, v1d.re, s1d.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, v1d.im, s1d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, v1d.re, s1d.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, v1d.im, s1d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(v1d.re, tempd.re);
    mpz_set(v1d.im, tempd.im);

    mpz_sub(v1d.re, v1d.re, l1d.re);
    mpz_sub(v1d.im, v1d.im, l1d.im);
    mpz_mul(temp1.re, u1d.re, l2d.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u1d.im, l2d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, u1d.re, l2d.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u1d.im, l2d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_add(v1d.re, v1d.re, temp1.re);
    mpz_add(v1d.im, v1d.im, temp1.im);
    mpz_mul(tempd.re, w3.re, v1d.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, w3.im, v1d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, w3.re, v1d.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, w3.im, v1d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(v1d.re, tempd.re);
    mpz_set(v1d.im, tempd.im);

    mpz_sub(v1d.re, v1d.re, m.v1.re);
    mpz_sub(v1d.im, v1d.im, m.v1.im);

    mpz_mul(v0d.re, u1d.re, u0d.re);
    mpz_mod(v0d.re, v0d.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u1d.im, u0d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(v0d.re, v0d.re, temp);
    mpz_mul(v0d.im, u1d.re, u0d.im);
    mpz_mod(v0d.im, v0d.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u1d.im, u0d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(v0d.im, v0d.im, temp);



    mpz_mul(tempd.re, v0d.re, s1d.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, v0d.im, s1d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, v0d.re, s1d.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, v0d.im, s1d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(v0d.re, tempd.re);
    mpz_set(v0d.im, tempd.im);

    mpz_add(v0d.re, v0d.re, l0d.re);
    mpz_add(v0d.im, v0d.im, l0d.im);
    mpz_mul(temp1.re, u0d.re, l2d.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u0d.im, l2d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, u0d.re, l2d.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u0d.im, l2d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_sub(v0d.re, temp1.re, v0d.re);
    mpz_sub(v0d.im, temp1.im, v0d.im);
    mpz_mul(tempd.re, v0d.re, w3.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, v0d.im, w3.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, v0d.re, w3.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, v0d.im, w3.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(v0d.re, tempd.re);
    mpz_set(v0d.im, tempd.im);

    mpz_sub(v0d.re, v0d.re, m.v0.re);
    mpz_sub(v0d.im, v0d.im, m.v0.im);

    ExtendedMumford ret(f, h, u1d, u0d, v1d, v0d);

    return ret;
}

ExtendedMumford ExtendedMumford::LangeDoubling() const{
    //std::cerr << "Lange Doubling." << std::endl;
    // 44M, 6S, I

    mpz_t temp;
    mpz_init(temp);
    ExtendedNumber temp1, temp2, tempd;

    // 1. v~ (=2v) と u の終結式を計算．
    // 3M, 2S
    ExtendedNumber v1t, v0t;
    mpz_add(v1t.re, v1.re, v1.re);
    mpz_add(v1t.im, v1.im, v1.im);
    mpz_add(v0t.re, v0.re, v0.re);
    mpz_add(v0t.im, v0.im, v0.im);

    ExtendedNumber w0, u1s, w1, w2, w3, r;
    mpz_mul(w0.re, v1.re, v1.re);
    mpz_mod(w0.re, w0.re, ExtendedNumber::CHARA);
    mpz_mul(temp, v1.im, v1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(w0.re, w0.re, temp);
    mpz_mul(w0.im, v1.re, v1.im);
    mpz_mod(w0.im, w0.im, ExtendedNumber::CHARA);
    mpz_mul(temp, v1.im, v1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(w0.im, w0.im, temp);



    mpz_mul(u1s.re, u1.re, u1.re);
    mpz_mod(u1s.re, u1s.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u1.im, u1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(u1s.re, u1s.re, temp);
    mpz_mul(u1s.im, u1.re, u1.im);
    mpz_mod(u1s.im, u1s.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u1.im, u1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(u1s.im, u1s.im, temp);



    mpz_add(w2.re, w0.re, w0.re);
    mpz_add(w2.im, w0.im, w0.im);
    mpz_add(w2.re, w2.re, w2.re);
    mpz_add(w2.im, w2.im, w2.im);
    mpz_mul(w3.re, u1.re, v1t.re);
    mpz_mod(w3.re, w3.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u1.im, v1t.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(w3.re, w3.re, temp);
    mpz_mul(w3.im, u1.re, v1t.im);
    mpz_mod(w3.im, w3.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u1.im, v1t.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(w3.im, w3.im, temp);



    mpz_mul(r.re, u0.re, w2.re);
    mpz_mod(r.re, r.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u0.im, w2.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(r.re, r.re, temp);
    mpz_mul(r.im, u0.re, w2.im);
    mpz_mod(r.im, r.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u0.im, w2.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(r.im, r.im, temp);



    mpz_sub(temp1.re, v0t.re, w3.re);
    mpz_sub(temp1.im, v0t.im, w3.im);
    mpz_mul(tempd.re, temp1.re, v0t.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, v0t.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, temp1.re, v0t.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, v0t.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(temp1.re, tempd.re);
    mpz_set(temp1.im, tempd.im);

    mpz_add(r.re, r.re, temp1.re);
    mpz_add(r.im, r.im, temp1.im);


    // 2. r の almost inverse を計算．
    //ExtendedNumber inv1d = -v1t;
    ExtendedNumber inv0d;
    mpz_sub(inv0d.re, v0t.re, w3.re);
    mpz_sub(inv0d.im, v0t.im, w3.im);
    //ExtendedNumber inv0d = v0t - w3;

    // 3. k を計算．
    // 11M
    //ExtendedNumber k4 = f6;
    ExtendedNumber k4u0, k4u1, k3, k3u0, k2, k1, k0, u1kd, k1d, k0d;
    mpz_mul(k4u0.re, this->f.coeff[6].re, u0.re);
    mpz_mod(k4u0.re, k4u0.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[6].im, u0.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(k4u0.re, k4u0.re, temp);
    mpz_mul(k4u0.im, this->f.coeff[6].re, u0.im);
    mpz_mod(k4u0.im, k4u0.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[6].im, u0.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(k4u0.im, k4u0.im, temp);



    mpz_mul(k4u1.re, this->f.coeff[6].re, u1.re);
    mpz_mod(k4u1.re, k4u1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[6].im, u1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(k4u1.re, k4u1.re, temp);
    mpz_mul(k4u1.im, this->f.coeff[6].re, u1.im);
    mpz_mod(k4u1.im, k4u1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[6].im, u1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(k4u1.im, k4u1.im, temp);



    mpz_sub(k3.re, this->f.coeff[5].re, k4u1.re);
    mpz_sub(k3.im, this->f.coeff[5].im, k4u1.im);
    mpz_mul(k3u0.re, k3.re, u0.re);
    mpz_mod(k3u0.re, k3u0.re, ExtendedNumber::CHARA);
    mpz_mul(temp, k3.im, u0.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(k3u0.re, k3u0.re, temp);
    mpz_mul(k3u0.im, k3.re, u0.im);
    mpz_mod(k3u0.im, k3u0.im, ExtendedNumber::CHARA);
    mpz_mul(temp, k3.im, u0.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(k3u0.im, k3u0.im, temp);



    mpz_mul(k2.re, k3.re, u1.re);
    mpz_mod(k2.re, k2.re, ExtendedNumber::CHARA);
    mpz_mul(temp, k3.im, u1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(k2.re, k2.re, temp);
    mpz_mul(k2.im, k3.re, u1.im);
    mpz_mod(k2.im, k2.im, ExtendedNumber::CHARA);
    mpz_mul(temp, k3.im, u1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(k2.im, k2.im, temp);



    mpz_sub(k2.re, this->f.coeff[4].re, k2.re);
    mpz_sub(k2.im, this->f.coeff[4].im, k2.im);
    mpz_sub(k2.re, k2.re, k4u0.re);
    mpz_sub(k2.im, k2.im, k4u0.im);
    mpz_mul(k1.re, k2.re, u1.re);
    mpz_mod(k1.re, k1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, k2.im, u1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(k1.re, k1.re, temp);
    mpz_mul(k1.im, k2.re, u1.im);
    mpz_mod(k1.im, k1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, k2.im, u1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(k1.im, k1.im, temp);



    mpz_sub(k1.re, this->f.coeff[3].re, k1.re);
    mpz_sub(k1.im, this->f.coeff[3].im, k1.im);
    mpz_sub(k1.re, k1.re, k3u0.re);
    mpz_sub(k1.im, k1.im, k3u0.im);
    mpz_mul(k0.re, k1.re, u1.re);
    mpz_mod(k0.re, k0.re, ExtendedNumber::CHARA);
    mpz_mul(temp, k1.im, u1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(k0.re, k0.re, temp);
    mpz_mul(k0.im, k1.re, u1.im);
    mpz_mod(k0.im, k0.im, ExtendedNumber::CHARA);
    mpz_mul(temp, k1.im, u1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(k0.im, k0.im, temp);


    mpz_mul(temp1.re, k2.re, u0.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, k2.im, u0.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, k2.re, u0.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, k2.im, u0.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_add(k0.re, k0.re, temp1.re);
    mpz_add(k0.im, k0.im, temp1.im);
    mpz_sub(k0.re, this->f.coeff[2].re, k0.re);
    mpz_sub(k0.im, this->f.coeff[2].im, k0.im);
    mpz_sub(k0.re, k0.re, w0.re);
    mpz_sub(k0.im, k0.im, w0.im);
    mpz_sub(u1kd.re, k3.re, k4u1.re);
    mpz_sub(u1kd.im, k3.im, k4u1.im);
    mpz_mul(tempd.re, u1.re, u1kd.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u1.im, u1kd.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, u1.re, u1kd.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u1.im, u1kd.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(u1kd.re, tempd.re);
    mpz_set(u1kd.im, tempd.im);


    mpz_add(k1d.re, k4u0.re, k4u0.re);
    mpz_add(k1d.im, k4u0.im, k4u0.im);
    mpz_sub(k1d.re, k1d.re, k2.re);
    mpz_sub(k1d.im, k1d.im, k2.im);
    mpz_mul(tempd.re, u1.re, k1d.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u1.im, k1d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, u1.re, k1d.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u1.im, k1d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(k1d.re, tempd.re);
    mpz_set(k1d.im, tempd.im);

    mpz_sub(temp1.re, k3.re, k4u1.re);
    mpz_sub(temp1.im, k3.im, k4u1.im);
    mpz_mul(tempd.re, temp1.re, u1s.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, u1s.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, temp1.re, u1s.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, u1s.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(temp1.re, tempd.re);
    mpz_set(temp1.im, tempd.im);

    mpz_add(k1d.re, k1d.re, temp1.re);
    mpz_add(k1d.im, k1d.im, temp1.im);
    mpz_sub(k1d.re, k1d.re, k3u0.re);
    mpz_sub(k1d.im, k1d.im, k3u0.im);
    mpz_add(k1d.re, k1d.re, k1.re);
    mpz_add(k1d.im, k1d.im, k1.im);

    mpz_add(k0d.re, u1kd.re, k4u0.re);
    mpz_add(k0d.im, u1kd.im, k4u0.im);
    mpz_sub(k0d.re, k0d.re, k2.re);
    mpz_sub(k0d.im, k0d.im, k2.im);
    mpz_mul(tempd.re, k0d.re, u0.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, k0d.im, u0.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, k0d.re, u0.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, k0d.im, u0.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(k0d.re, tempd.re);
    mpz_set(k0d.im, tempd.im);

    mpz_add(k0d.re, k0d.re, k0.re);
    mpz_add(k0d.im, k0d.im, k0.im);

    // 4. s' を計算．
    // 5M
    mpz_mul(w0.re, k0d.re, inv0d.re);
    mpz_mod(w0.re, w0.re, ExtendedNumber::CHARA);
    mpz_mul(temp, k0d.im, inv0d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(w0.re, w0.re, temp);
    mpz_mul(w0.im, k0d.re, inv0d.im);
    mpz_mod(w0.im, w0.im, ExtendedNumber::CHARA);
    mpz_mul(temp, k0d.im, inv0d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(w0.im, w0.im, temp);



    mpz_mul(w1.re, k1d.re, v1t.re);
    mpz_mod(w1.re, w1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, k1d.im, v1t.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(w1.re, w1.re, temp);
    mpz_mul(w1.im, k1d.re, v1t.im);
    mpz_mod(w1.im, w1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, k1d.im, v1t.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(w1.im, w1.im, temp);



    mpz_neg(w1.re, w1.re);
    mpz_neg(w1.im, w1.im);
    ExtendedNumber s0d, s1d;
    mpz_sub(s1d.re, inv0d.re, v1t.re);
    mpz_sub(s1d.im, inv0d.im, v1t.im);
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

    mpz_sub(s1d.re, s1d.re, w0.re);
    mpz_sub(s1d.im, s1d.im, w0.im);
    mpz_add(temp1.re, ExtendedNumber::ONE().re, u1.re);
    mpz_add(temp1.im, ExtendedNumber::ONE().im, u1.im);
    mpz_mul(tempd.re, temp1.re, w1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, w1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, temp1.re, w1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, temp1.im, w1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(temp1.re, tempd.re);
    mpz_set(temp1.im, tempd.im);

    mpz_sub(s1d.re, s1d.re, temp1.re);
    mpz_sub(s1d.im, s1d.im, temp1.im);
    mpz_mul(s0d.re, u0.re, w1.re);
    mpz_mod(s0d.re, s0d.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u0.im, w1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(s0d.re, s0d.re, temp);
    mpz_mul(s0d.im, u0.re, w1.im);
    mpz_mod(s0d.im, s0d.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u0.im, w1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(s0d.im, s0d.im, temp);



    mpz_sub(s0d.re, w0.re, s0d.re);
    mpz_sub(s0d.im, w0.im, s0d.im);

    if(s1d.isZero()){
        std::cerr << "Special case." << std::endl;
        // サブルーチン
        ExtendedMumford ret(f, h);
        return ret;
    }

    // 5. u' を計算．
    // 8M, 3S
    ExtendedNumber rs, u2d, u1d, u0d;
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



    mpz_mul(u2d.re, s1d.re, s1d.re);
    mpz_mod(u2d.re, u2d.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, s1d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(u2d.re, u2d.re, temp);
    mpz_mul(u2d.im, s1d.re, s1d.im);
    mpz_mod(u2d.im, u2d.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, s1d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(u2d.im, u2d.im, temp);



    mpz_mul(temp1.re, rs.re, this->f.coeff[6].re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, rs.im, this->f.coeff[6].im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, rs.re, this->f.coeff[6].im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, rs.im, this->f.coeff[6].re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_sub(u2d.re, u2d.re, temp1.re);
    mpz_sub(u2d.im, u2d.im, temp1.im);
    mpz_mul(u1d.re, s1d.re, s0d.re);
    mpz_mod(u1d.re, u1d.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, s0d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(u1d.re, u1d.re, temp);
    mpz_mul(u1d.im, s1d.re, s0d.im);
    mpz_mod(u1d.im, u1d.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, s0d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(u1d.im, u1d.im, temp);



    mpz_add(u1d.re, u1d.re, u1d.re);
    mpz_add(u1d.im, u1d.im, u1d.im);
    mpz_sub(temp1.re, this->f.coeff[5].re, k4u1.re);
    mpz_sub(temp1.im, this->f.coeff[5].im, k4u1.im);
    mpz_sub(temp1.re, temp1.re, k4u1.re);
    mpz_sub(temp1.im, temp1.im, k4u1.im);
    mpz_mul(tempd.re, rs.re, temp1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, rs.im, temp1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, rs.re, temp1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, rs.im, temp1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(temp1.re, tempd.re);
    mpz_set(temp1.im, tempd.im);

    mpz_sub(u1d.re, u1d.re, temp1.re);
    mpz_sub(u1d.im, u1d.im, temp1.im);

    mpz_mul(u0d.re, v1.re, s1d.re);
    mpz_mod(u0d.re, u0d.re, ExtendedNumber::CHARA);
    mpz_mul(temp, v1.im, s1d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(u0d.re, u0d.re, temp);
    mpz_mul(u0d.im, v1.re, s1d.im);
    mpz_mod(u0d.im, u0d.im, ExtendedNumber::CHARA);
    mpz_mul(temp, v1.im, s1d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(u0d.im, u0d.im, temp);



    mpz_mul(tempd.re, u0d.re, r.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u0d.im, r.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, u0d.re, r.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u0d.im, r.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(u0d.re, tempd.re);
    mpz_set(u0d.im, tempd.im);

    mpz_add(u0d.re, u0d.re, u0d.re);
    mpz_add(u0d.im, u0d.im, u0d.im);
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



    mpz_add(u0d.re, u0d.re, temp1.re);
    mpz_add(u0d.im, u0d.im, temp1.im);

    mpz_sub(temp1.re, this->f.coeff[5].re, k4u1.re);
    mpz_sub(temp1.im, this->f.coeff[5].im, k4u1.im);
    mpz_sub(temp1.re, temp1.re, k4u1.re);
    mpz_sub(temp1.im, temp1.im, k4u1.im);
    mpz_add(temp1.re, temp1.re, temp1.re);
    mpz_add(temp1.im, temp1.im, temp1.im);
    mpz_mul(tempd.re, u1.re, temp1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u1.im, temp1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, u1.re, temp1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u1.im, temp1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(temp1.re, tempd.re);
    mpz_set(temp1.im, tempd.im);

    mpz_add(temp2.re, u0.re, u0.re);
    mpz_add(temp2.im, u0.im, u0.im);
    mpz_add(temp2.re, temp2.re, u1s.re);
    mpz_add(temp2.im, temp2.im, u1s.im);
    mpz_mul(tempd.re, temp2.re, this->f.coeff[6].re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, temp2.im, this->f.coeff[6].im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, temp2.re, this->f.coeff[6].im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, temp2.im, this->f.coeff[6].re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(temp2.re, tempd.re);
    mpz_set(temp2.im, tempd.im);

    mpz_sub(temp2.re, this->f.coeff[4].re, temp2.re);
    mpz_sub(temp2.im, this->f.coeff[4].im, temp2.im);
    mpz_sub(temp2.re, temp2.re, temp1.re);
    mpz_sub(temp2.im, temp2.im, temp1.im);
    mpz_mul(tempd.re, temp2.re, rs.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, temp2.im, rs.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, temp2.re, rs.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, temp2.im, rs.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(temp2.re, tempd.re);
    mpz_set(temp2.im, tempd.im);

    mpz_sub(u0d.re, u0d.re, temp2.re);
    mpz_sub(u0d.im, u0d.im, temp2.im);

    // 6. r と u2d の逆元を計算．
    // 3M, I
    mpz_mul(w0.re, r.re, u2d.re);
    mpz_mod(w0.re, w0.re, ExtendedNumber::CHARA);
    mpz_mul(temp, r.im, u2d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(w0.re, w0.re, temp);
    mpz_mul(w0.im, r.re, u2d.im);
    mpz_mod(w0.im, w0.im, ExtendedNumber::CHARA);
    mpz_mul(temp, r.im, u2d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(w0.im, w0.im, temp);



    mpz_t denom;
    mpz_init(denom);
    mpz_mul(temp, w0.re, w0.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_mul(denom, w0.im, w0.im);
    mpz_mod(denom, denom, ExtendedNumber::CHARA);
    mpz_add(denom, denom, temp);
    constant_invert(denom, denom, ExtendedNumber::CHARA);
    mpz_mul(w0.re, w0.re, denom);
    mpz_mod(w0.re, w0.re, ExtendedNumber::CHARA);
    mpz_mul(w0.im, w0.im, denom);
    mpz_mod(w0.im, w0.im, ExtendedNumber::CHARA);
    mpz_neg(w0.im, w0.im);

    mpz_mul(w1.re, w0.re, r.re);
    mpz_mod(w1.re, w1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, w0.im, r.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(w1.re, w1.re, temp);
    mpz_mul(w1.im, w0.re, r.im);
    mpz_mod(w1.im, w1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, w0.im, r.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(w1.im, w1.im, temp);



    mpz_mul(w2.re, w0.re, u2d.re);
    mpz_mod(w2.re, w2.re, ExtendedNumber::CHARA);
    mpz_mul(temp, w0.im, u2d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(w2.re, w2.re, temp);
    mpz_mul(w2.im, w0.re, u2d.im);
    mpz_mod(w2.im, w2.im, ExtendedNumber::CHARA);
    mpz_mul(temp, w0.im, u2d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(w2.im, w2.im, temp);




    // 7. u' を計算．
    // 2M
    mpz_mul(tempd.re, u1d.re, w1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u1d.im, w1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, u1d.re, w1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u1d.im, w1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(u1d.re, tempd.re);
    mpz_set(u1d.im, tempd.im);

    mpz_mul(tempd.re, u0d.re, w1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u0d.im, w1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, u0d.re, w1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u0d.im, w1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(u0d.re, tempd.re);
    mpz_set(u0d.im, tempd.im);


    // 8. v' を計算．
    // 12M, 1S
    ExtendedNumber l3, l2, l1, l0, v1d, v0d;
    mpz_mul(l3.re, s1d.re, w2.re);
    mpz_mod(l3.re, l3.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, w2.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(l3.re, l3.re, temp);
    mpz_mul(l3.im, s1d.re, w2.im);
    mpz_mod(l3.im, l3.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, w2.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(l3.im, l3.im, temp);



    mpz_mul(l2.re, s1d.re, u1.re);
    mpz_mod(l2.re, l2.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, u1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(l2.re, l2.re, temp);
    mpz_mul(l2.im, s1d.re, u1.im);
    mpz_mod(l2.im, l2.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, u1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(l2.im, l2.im, temp);



    mpz_add(l2.re, l2.re, s0d.re);
    mpz_add(l2.im, l2.im, s0d.im);
    mpz_mul(tempd.re, l2.re, w2.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, w2.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, l2.re, w2.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, w2.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(l2.re, tempd.re);
    mpz_set(l2.im, tempd.im);

    mpz_mul(l1.re, s1d.re, u0.re);
    mpz_mod(l1.re, l1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, u0.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(l1.re, l1.re, temp);
    mpz_mul(l1.im, s1d.re, u0.im);
    mpz_mod(l1.im, l1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s1d.im, u0.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(l1.im, l1.im, temp);



    mpz_mul(temp1.re, s0d.re, u1.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s0d.im, u1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, s0d.re, u1.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s0d.im, u1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_add(l1.re, l1.re, temp1.re);
    mpz_add(l1.im, l1.im, temp1.im);
    mpz_mul(tempd.re, l1.re, w2.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l1.im, w2.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, l1.re, w2.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l1.im, w2.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(l1.re, tempd.re);
    mpz_set(l1.im, tempd.im);

    mpz_mul(l0.re, s0d.re, u0.re);
    mpz_mod(l0.re, l0.re, ExtendedNumber::CHARA);
    mpz_mul(temp, s0d.im, u0.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(l0.re, l0.re, temp);
    mpz_mul(l0.im, s0d.re, u0.im);
    mpz_mod(l0.im, l0.im, ExtendedNumber::CHARA);
    mpz_mul(temp, s0d.im, u0.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(l0.im, l0.im, temp);



    mpz_mul(tempd.re, l0.re, w2.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l0.im, w2.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, l0.re, w2.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l0.im, w2.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(l0.re, tempd.re);
    mpz_set(l0.im, tempd.im);


    mpz_mul(v1d.re, u1d.re, u1d.re);
    mpz_mod(v1d.re, v1d.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u1d.im, u1d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(v1d.re, v1d.re, temp);
    mpz_mul(v1d.im, u1d.re, u1d.im);
    mpz_mod(v1d.im, v1d.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u1d.im, u1d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(v1d.im, v1d.im, temp);



    mpz_sub(v1d.re, v1d.re, u0d.re);
    mpz_sub(v1d.im, v1d.im, u0d.im);
    mpz_mul(tempd.re, v1d.re, l3.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, v1d.im, l3.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, v1d.re, l3.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, v1d.im, l3.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(v1d.re, tempd.re);
    mpz_set(v1d.im, tempd.im);

    mpz_mul(temp1.re, l2.re, u1d.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, u1d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, l2.re, u1d.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, u1d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_sub(v1d.re, v1d.re, temp1.re);
    mpz_sub(v1d.im, v1d.im, temp1.im);
    mpz_add(v1d.re, v1d.re, l1.re);
    mpz_add(v1d.im, v1d.im, l1.im);
    mpz_add(v1d.re, v1d.re, v1.re);
    mpz_add(v1d.im, v1d.im, v1.im);

    mpz_mul(v0d.re, l3.re, u1d.re);
    mpz_mod(v0d.re, v0d.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l3.im, u1d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(v0d.re, v0d.re, temp);
    mpz_mul(v0d.im, l3.re, u1d.im);
    mpz_mod(v0d.im, v0d.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l3.im, u1d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(v0d.im, v0d.im, temp);



    mpz_sub(v0d.re, v0d.re, l2.re);
    mpz_sub(v0d.im, v0d.im, l2.im);
    mpz_mul(tempd.re, v0d.re, u0d.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, v0d.im, u0d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, v0d.re, u0d.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, v0d.im, u0d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(v0d.re, tempd.re);
    mpz_set(v0d.im, tempd.im);

    mpz_add(v0d.re, v0d.re, l0.re);
    mpz_add(v0d.im, v0d.im, l0.im);
    mpz_add(v0d.re, v0d.re, v0.re);
    mpz_add(v0d.im, v0d.im, v0.im);
    mpz_neg(v1d.re, v1d.re);
    mpz_neg(v1d.im, v1d.im);
    mpz_neg(v0d.re, v0d.re);
    mpz_neg(v0d.im, v0d.im);

    ExtendedMumford ret(f, h, u1d, u0d, v1d, v0d);
    return ret;
}

ExtendedMumford ExtendedMumford::CostelloDoubling() const{
    //std::cerr << "Costello Doubling." << std::endl;
    // 32M, 6S, I

    mpz_t temp;
    mpz_init(temp);
    ExtendedNumber temp1, temp2, temp3, temp4, temp5, tempd;

    ExtendedNumber vv, va;
    mpz_mul(vv.re, this->v1.re, this->v1.re);
    mpz_mod(vv.re, vv.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->v1.im, this->v1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(vv.re, vv.re, temp);
    mpz_mul(vv.im, this->v1.re, this->v1.im);
    mpz_mod(vv.im, vv.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->v1.im, this->v1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(vv.im, vv.im, temp);



    mpz_add(va.re, this->v1.re, this->u1.re);
    mpz_add(va.im, this->v1.im, this->u1.im);
    mpz_mul(tempd.re, va.re, va.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, va.im, va.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, va.re, va.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, va.im, va.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(va.re, tempd.re);
    mpz_set(va.im, tempd.im);

    mpz_sub(va.re, va.re, vv.re);
    mpz_sub(va.im, va.im, vv.im);
    mpz_sub(va.re, va.re, this->U1.re);
    mpz_sub(va.im, va.im, this->U1.im);

    ExtendedNumber M1, M2, M3, M4;
    mpz_sub(M1.re, this->v0.re, va.re);
    mpz_sub(M1.im, this->v0.im, va.im);
    mpz_add(M1.re, M1.re, M1.re);
    mpz_add(M1.im, M1.im, M1.im);
    mpz_add(M2.re, this->U1.re, this->U1.re);
    mpz_add(M2.im, this->U1.im, this->U1.im);
    mpz_add(M2.re, M2.re, this->u0.re);
    mpz_add(M2.im, M2.im, this->u0.im);
    mpz_mul(tempd.re, M2.re, this->v1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, M2.im, this->v1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, M2.re, this->v1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, M2.im, this->v1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(M2.re, tempd.re);
    mpz_set(M2.im, tempd.im);

    mpz_add(M2.re, M2.re, M2.re);
    mpz_add(M2.im, M2.im, M2.im);
    mpz_add(M3.re, this->v1.re, this->v1.re);
    mpz_add(M3.im, this->v1.im, this->v1.im);
    mpz_neg(M3.re, M3.re);
    mpz_neg(M3.im, M3.im);
    mpz_add(M4.re, va.re, this->v0.re);
    mpz_add(M4.im, va.im, this->v0.im);
    mpz_add(M4.re, M4.re, this->v0.re);
    mpz_add(M4.im, M4.im, this->v0.im);

    mpz_mul(temp2.re, this->f.coeff[6].re, this->u0.re);
    mpz_mod(temp2.re, temp2.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[6].im, this->u0.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp2.re, temp2.re, temp);
    mpz_mul(temp2.im, this->f.coeff[6].re, this->u0.im);
    mpz_mod(temp2.im, temp2.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[6].im, this->u0.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp2.im, temp2.im, temp);



    mpz_mul(temp3.re, this->f.coeff[6].re, this->U1.re);
    mpz_mod(temp3.re, temp3.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[6].im, this->U1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp3.re, temp3.re, temp);
    mpz_mul(temp3.im, this->f.coeff[6].re, this->U1.im);
    mpz_mod(temp3.im, temp3.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[6].im, this->U1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp3.im, temp3.im, temp);



    mpz_mul(temp4.re, this->f.coeff[5].re, this->u0.re);
    mpz_mod(temp4.re, temp4.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[5].im, this->u0.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp4.re, temp4.re, temp);
    mpz_mul(temp4.im, this->f.coeff[5].re, this->u0.im);
    mpz_mod(temp4.im, temp4.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[5].im, this->u0.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp4.im, temp4.im, temp);



    mpz_mul(temp5.re, this->f.coeff[5].re, this->u1.re);
    mpz_mod(temp5.re, temp5.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[5].im, this->u1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp5.re, temp5.re, temp);
    mpz_mul(temp5.im, this->f.coeff[5].re, this->u1.im);
    mpz_mod(temp5.im, temp5.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[5].im, this->u1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp5.im, temp5.im, temp);




    ExtendedNumber z11, z12;
    mpz_sub(z11.re, temp5.re, temp3.re);
    mpz_sub(z11.im, temp5.im, temp3.im);
    mpz_add(z11.re, z11.re, z11.re);
    mpz_add(z11.im, z11.im, z11.im);
    mpz_sub(z11.re, z11.re, temp3.re);
    mpz_sub(z11.im, z11.im, temp3.im);
    mpz_sub(z11.re, z11.re, this->f.coeff[4].re);
    mpz_sub(z11.im, z11.im, this->f.coeff[4].im);
    mpz_mul(tempd.re, z11.re, this->U1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, z11.im, this->U1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, z11.re, this->U1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, z11.im, this->U1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(z11.re, tempd.re);
    mpz_set(z11.im, tempd.im);


    mpz_add(z12.re, temp5.re, temp2.re);
    mpz_add(z12.im, temp5.im, temp2.im);
    mpz_add(z12.re, z12.re, z12.re);
    mpz_add(z12.im, z12.im, z12.im);
    mpz_add(z12.re, z12.re, temp2.re);
    mpz_add(z12.im, z12.im, temp2.im);
    mpz_sub(z12.re, z12.re, this->f.coeff[4].re);
    mpz_sub(z12.im, z12.im, this->f.coeff[4].im);
    mpz_sub(z12.re, z12.re, this->f.coeff[4].re);
    mpz_sub(z12.im, z12.im, this->f.coeff[4].im);
    mpz_mul(tempd.re, z12.re, this->u0.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, z12.im, this->u0.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, z12.re, this->u0.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, z12.im, this->u0.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(z12.re, tempd.re);
    mpz_set(z12.im, tempd.im);


    ExtendedNumber z1, z2;
    mpz_add(z1.re, z11.re, z12.re);
    mpz_add(z1.im, z11.im, z12.im);
    mpz_sub(z1.re, z1.re, vv.re);
    mpz_sub(z1.im, z1.im, vv.im);
    mpz_add(z1.re, z1.re, this->f.coeff[2].re);
    mpz_add(z1.im, z1.im, this->f.coeff[2].im);
    mpz_sub(z2.re, temp2.re, temp3.re);
    mpz_sub(z2.im, temp2.im, temp3.im);
    mpz_add(z2.re, z2.re, temp5.re);
    mpz_add(z2.im, z2.im, temp5.im);
    mpz_sub(z2.re, z2.re, this->f.coeff[4].re);
    mpz_sub(z2.im, z2.im, this->f.coeff[4].im);
    mpz_add(z2.re, z2.re, z2.re);
    mpz_add(z2.im, z2.im, z2.im);
    mpz_add(z2.re, z2.re, temp2.re);
    mpz_add(z2.im, z2.im, temp2.im);
    mpz_add(z2.re, z2.re, temp2.re);
    mpz_add(z2.im, z2.im, temp2.im);
    mpz_add(z2.re, z2.re, temp2.re);
    mpz_add(z2.im, z2.im, temp2.im);
    mpz_add(z2.re, z2.re, temp2.re);
    mpz_add(z2.im, z2.im, temp2.im);
    mpz_sub(z2.re, z2.re, temp3.re);
    mpz_sub(z2.im, z2.im, temp3.im);
    mpz_sub(z2.re, z2.re, temp3.re);
    mpz_sub(z2.im, z2.im, temp3.im);
    mpz_add(z2.re, z2.re, temp5.re);
    mpz_add(z2.im, z2.im, temp5.im);
    mpz_mul(tempd.re, z2.re, this->u1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, z2.im, this->u1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, z2.re, this->u1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, z2.im, this->u1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(z2.re, tempd.re);
    mpz_set(z2.im, tempd.im);

    mpz_sub(z2.re, z2.re, temp4.re);
    mpz_sub(z2.im, z2.im, temp4.im);
    mpz_sub(z2.re, z2.re, temp4.re);
    mpz_sub(z2.im, z2.im, temp4.im);
    mpz_add(z2.re, z2.re, this->f.coeff[3].re);
    mpz_add(z2.im, z2.im, this->f.coeff[3].im);

    ExtendedNumber t1, t2, t3, t4;
    mpz_sub(t1.re, M2.re, z1.re);
    mpz_sub(t1.im, M2.im, z1.im);
    mpz_sub(temp1.re, z2.re, M1.re);
    mpz_sub(temp1.im, z2.im, M1.im);
    mpz_mul(tempd.re, t1.re, temp1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, t1.im, temp1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, t1.re, temp1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, t1.im, temp1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(t1.re, tempd.re);
    mpz_set(t1.im, tempd.im);


    mpz_add(t2.re, z1.re, M2.re);
    mpz_add(t2.im, z1.im, M2.im);
    mpz_add(temp1.re, z2.re, M1.re);
    mpz_add(temp1.im, z2.im, M1.im);
    mpz_mul(tempd.re, t2.re, temp1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, t2.im, temp1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, t2.re, temp1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, t2.im, temp1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(t2.re, tempd.re);
    mpz_set(t2.im, tempd.im);

    mpz_neg(t2.re, t2.re);
    mpz_neg(t2.im, t2.im);

    mpz_sub(t3.re, M4.re, z1.re);
    mpz_sub(t3.im, M4.im, z1.im);
    mpz_sub(temp1.re, z2.re, M3.re);
    mpz_sub(temp1.im, z2.im, M3.im);
    mpz_mul(tempd.re, t3.re, temp1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, t3.im, temp1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, t3.re, temp1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, t3.im, temp1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(t3.re, tempd.re);
    mpz_set(t3.im, tempd.im);


    mpz_add(t4.re, M4.re, z1.re);
    mpz_add(t4.im, M4.im, z1.im);
    mpz_add(temp1.re, z2.re, M3.re);
    mpz_add(temp1.im, z2.im, M3.im);
    mpz_mul(tempd.re, t4.re, temp1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, t4.im, temp1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, t4.re, temp1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, t4.im, temp1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(t4.re, tempd.re);
    mpz_set(t4.im, tempd.im);

    mpz_neg(t4.re, t4.re);
    mpz_neg(t4.im, t4.im);

    ExtendedNumber l2_num, l3_num;
    mpz_sub(l2_num.re, t1.re, t2.re);
    mpz_sub(l2_num.im, t1.im, t2.im);
    mpz_sub(l3_num.re, t3.re, t4.re);
    mpz_sub(l3_num.im, t3.im, t4.im);

    ExtendedNumber d;
    mpz_sub(d.re, M4.re, M2.re);
    mpz_sub(d.im, M4.im, M2.im);
    mpz_add(temp1.re, M1.re, M3.re);
    mpz_add(temp1.im, M1.im, M3.im);
    mpz_mul(tempd.re, d.re, temp1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, d.im, temp1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, d.re, temp1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, d.im, temp1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(d.re, tempd.re);
    mpz_set(d.im, tempd.im);

    mpz_add(d.re, d.re, d.re);
    mpz_add(d.im, d.im, d.im);
    mpz_add(d.re, d.re, t3.re);
    mpz_add(d.im, d.im, t3.im);
    mpz_add(d.re, d.re, t4.re);
    mpz_add(d.im, d.im, t4.im);
    mpz_sub(d.re, d.re, t1.re);
    mpz_sub(d.im, d.im, t1.im);
    mpz_sub(d.re, d.re, t2.re);
    mpz_sub(d.im, d.im, t2.im);

    ExtendedNumber A, B, C, d_inv, d_shifted_inv;
    mpz_mul(A.re, d.re, d.re);
    mpz_mod(A.re, A.re, ExtendedNumber::CHARA);
    mpz_mul(temp, d.im, d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(A.re, A.re, temp);
    mpz_mul(A.im, d.re, d.im);
    mpz_mod(A.im, A.im, ExtendedNumber::CHARA);
    mpz_mul(temp, d.im, d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(A.im, A.im, temp);



    mpz_mul(temp1.re, this->f.coeff[6].re, A.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[6].im, A.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, this->f.coeff[6].re, A.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, this->f.coeff[6].im, A.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_mul(B.re, l3_num.re, l3_num.re);
    mpz_mod(B.re, B.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l3_num.im, l3_num.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(B.re, B.re, temp);
    mpz_mul(B.im, l3_num.re, l3_num.im);
    mpz_mod(B.im, B.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l3_num.im, l3_num.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(B.im, B.im, temp);



    mpz_sub(B.re, B.re, temp1.re);
    mpz_sub(B.im, B.im, temp1.im);
    mpz_mul(C.re, d.re, B.re);
    mpz_mod(C.re, C.re, ExtendedNumber::CHARA);
    mpz_mul(temp, d.im, B.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(C.re, C.re, temp);
    mpz_mul(C.im, d.re, B.im);
    mpz_mod(C.im, C.im, ExtendedNumber::CHARA);
    mpz_mul(temp, d.im, B.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(C.im, C.im, temp);



    mpz_t denom;
    mpz_init(denom);
    mpz_mul(temp, C.re, C.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_mul(denom, C.im, C.im);
    mpz_mod(denom, denom, ExtendedNumber::CHARA);
    mpz_add(denom, denom, temp);
    constant_invert(denom, denom, ExtendedNumber::CHARA);
    mpz_mul(C.re, C.re, denom);
    mpz_mod(C.re, C.re, ExtendedNumber::CHARA);
    mpz_mul(C.im, C.im, denom);
    mpz_mod(C.im, C.im, ExtendedNumber::CHARA);
    mpz_neg(C.im, C.im);

    mpz_mul(d_inv.re, B.re, C.re);
    mpz_mod(d_inv.re, d_inv.re, ExtendedNumber::CHARA);
    mpz_mul(temp, B.im, C.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(d_inv.re, d_inv.re, temp);
    mpz_mul(d_inv.im, B.re, C.im);
    mpz_mod(d_inv.im, d_inv.im, ExtendedNumber::CHARA);
    mpz_mul(temp, B.im, C.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(d_inv.im, d_inv.im, temp);



    mpz_mul(d_shifted_inv.re, d.re, A.re);
    mpz_mod(d_shifted_inv.re, d_shifted_inv.re, ExtendedNumber::CHARA);
    mpz_mul(temp, d.im, A.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(d_shifted_inv.re, d_shifted_inv.re, temp);
    mpz_mul(d_shifted_inv.im, d.re, A.im);
    mpz_mod(d_shifted_inv.im, d_shifted_inv.im, ExtendedNumber::CHARA);
    mpz_mul(temp, d.im, A.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(d_shifted_inv.im, d_shifted_inv.im, temp);



    mpz_mul(tempd.re, d_shifted_inv.re, C.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, d_shifted_inv.im, C.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, d_shifted_inv.re, C.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, d_shifted_inv.im, C.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(d_shifted_inv.re, tempd.re);
    mpz_set(d_shifted_inv.im, tempd.im);


    ExtendedNumber l2, l3;
    mpz_mul(l2.re, l2_num.re, d_inv.re);
    mpz_mod(l2.re, l2.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l2_num.im, d_inv.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(l2.re, l2.re, temp);
    mpz_mul(l2.im, l2_num.re, d_inv.im);
    mpz_mod(l2.im, l2.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l2_num.im, d_inv.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(l2.im, l2.im, temp);



    mpz_mul(l3.re, l3_num.re, d_inv.re);
    mpz_mod(l3.re, l3.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l3_num.im, d_inv.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(l3.re, l3.re, temp);
    mpz_mul(l3.im, l3_num.re, d_inv.im);
    mpz_mod(l3.im, l3.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l3_num.im, d_inv.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(l3.im, l3.im, temp);




    mpz_mul(temp1.re, l2.re, l3.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, l3.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, l2.re, l3.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, l3.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);




    ExtendedNumber u1d;
    mpz_sub(u1d.re, this->f.coeff[5].re, temp1.re);
    mpz_sub(u1d.im, this->f.coeff[5].im, temp1.im);
    mpz_sub(u1d.re, u1d.re, temp1.re);
    mpz_sub(u1d.im, u1d.im, temp1.im);
    mpz_mul(tempd.re, u1d.re, d_shifted_inv.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u1d.im, d_shifted_inv.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, u1d.re, d_shifted_inv.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u1d.im, d_shifted_inv.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(u1d.re, tempd.re);
    mpz_set(u1d.im, tempd.im);

    mpz_add(u1d.re, u1d.re, this->u1.re);
    mpz_add(u1d.im, u1d.im, this->u1.im);
    mpz_add(u1d.re, u1d.re, this->u1.re);
    mpz_add(u1d.im, u1d.im, this->u1.im);
    mpz_neg(u1d.re, u1d.re);
    mpz_neg(u1d.im, u1d.im);

    ExtendedNumber u0d;
    mpz_sub(u0d.re, this->u0.re, this->U1.re);
    mpz_sub(u0d.im, this->u0.im, this->U1.im);
    mpz_mul(tempd.re, u0d.re, l3.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u0d.im, l3.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, u0d.re, l3.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u0d.im, l3.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(u0d.re, tempd.re);
    mpz_set(u0d.im, tempd.im);

    mpz_mul(temp1.re, l2.re, this->u1.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, this->u1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, l2.re, this->u1.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, this->u1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_add(u0d.re, u0d.re, temp1.re);
    mpz_add(u0d.im, u0d.im, temp1.im);
    mpz_add(u0d.re, u0d.re, this->v1.re);
    mpz_add(u0d.im, u0d.im, this->v1.im);
    mpz_mul(tempd.re, u0d.re, l3.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u0d.im, l3.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, u0d.re, l3.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u0d.im, l3.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(u0d.re, tempd.re);
    mpz_set(u0d.im, tempd.im);

    mpz_add(u0d.re, u0d.re, u0d.re);
    mpz_add(u0d.im, u0d.im, u0d.im);
    mpz_mul(temp1.re, l2.re, l2.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, l2.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, l2.re, l2.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, l2.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_add(u0d.re, u0d.re, temp1.re);
    mpz_add(u0d.im, u0d.im, temp1.im);
    mpz_sub(u0d.re, u0d.re, this->f.coeff[4].re);
    mpz_sub(u0d.im, u0d.im, this->f.coeff[4].im);
    mpz_mul(tempd.re, u0d.re, d_shifted_inv.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u0d.im, d_shifted_inv.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, u0d.re, d_shifted_inv.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u0d.im, d_shifted_inv.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(u0d.re, tempd.re);
    mpz_set(u0d.im, tempd.im);

    mpz_mul(temp1.re, u1d.re, this->u1.re);
    mpz_mod(temp1.re, temp1.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u1d.im, this->u1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(temp1.re, temp1.re, temp);
    mpz_mul(temp1.im, u1d.re, this->u1.im);
    mpz_mod(temp1.im, temp1.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u1d.im, this->u1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(temp1.im, temp1.im, temp);



    mpz_sub(u0d.re, u0d.re, temp1.re);
    mpz_sub(u0d.im, u0d.im, temp1.im);
    mpz_sub(u0d.re, u0d.re, temp1.re);
    mpz_sub(u0d.im, u0d.im, temp1.im);
    mpz_sub(u0d.re, u0d.re, this->u0.re);
    mpz_sub(u0d.im, u0d.im, this->u0.im);
    mpz_sub(u0d.re, u0d.re, this->u0.re);
    mpz_sub(u0d.im, u0d.im, this->u0.im);
    mpz_sub(u0d.re, u0d.re, this->U1.re);
    mpz_sub(u0d.im, u0d.im, this->U1.im);

    ExtendedNumber U1d, U0d;
    mpz_mul(U1d.re, u1d.re, u1d.re);
    mpz_mod(U1d.re, U1d.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u1d.im, u1d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(U1d.re, U1d.re, temp);
    mpz_mul(U1d.im, u1d.re, u1d.im);
    mpz_mod(U1d.im, U1d.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u1d.im, u1d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(U1d.im, U1d.im, temp);



    mpz_mul(U0d.re, u1d.re, u0d.re);
    mpz_mod(U0d.re, U0d.re, ExtendedNumber::CHARA);
    mpz_mul(temp, u1d.im, u0d.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(U0d.re, U0d.re, temp);
    mpz_mul(U0d.im, u1d.re, u0d.im);
    mpz_mod(U0d.im, U0d.im, ExtendedNumber::CHARA);
    mpz_mul(temp, u1d.im, u0d.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(U0d.im, U0d.im, temp);




    ExtendedNumber v1d, v0d;
    mpz_sub(temp1.re, u1d.re, this->u1.re);
    mpz_sub(temp1.im, u1d.im, this->u1.im);
    mpz_mul(tempd.re, l2.re, temp1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, temp1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, l2.re, temp1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, temp1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(temp1.re, tempd.re);
    mpz_set(temp1.im, tempd.im);


    mpz_sub(v1d.re, u0d.re, U1d.re);
    mpz_sub(v1d.im, u0d.im, U1d.im);
    mpz_add(v1d.re, v1d.re, this->U1.re);
    mpz_add(v1d.im, v1d.im, this->U1.im);
    mpz_sub(v1d.re, v1d.re, this->u0.re);
    mpz_sub(v1d.im, v1d.im, this->u0.im);
    mpz_mul(tempd.re, v1d.re, l3.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, v1d.im, l3.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, v1d.re, l3.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, v1d.im, l3.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(v1d.re, tempd.re);
    mpz_set(v1d.im, tempd.im);

    mpz_add(v1d.re, v1d.re, temp1.re);
    mpz_add(v1d.im, v1d.im, temp1.im);
    mpz_sub(v1d.re, v1d.re, this->v1.re);
    mpz_sub(v1d.im, v1d.im, this->v1.im);

    mpz_sub(temp1.re, u0d.re, this->u0.re);
    mpz_sub(temp1.im, u0d.im, this->u0.im);
    mpz_mul(tempd.re, l2.re, temp1.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, temp1.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, l2.re, temp1.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, l2.im, temp1.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(temp1.re, tempd.re);
    mpz_set(temp1.im, tempd.im);


    mpz_sub(v0d.re, this->U0.re, U0d.re);
    mpz_sub(v0d.im, this->U0.im, U0d.im);
    mpz_mul(tempd.re, v0d.re, l3.re);
    mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);
    mpz_mul(temp, v0d.im, l3.im);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_sub(tempd.re, tempd.re, temp);
    mpz_mul(tempd.im, v0d.re, l3.im);
    mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);
    mpz_mul(temp, v0d.im, l3.re);
    mpz_mod(temp, temp, ExtendedNumber::CHARA);
    mpz_add(tempd.im, tempd.im, temp);
    mpz_set(v0d.re, tempd.re);
    mpz_set(v0d.im, tempd.im);

    mpz_add(v0d.re, v0d.re, temp1.re);
    mpz_add(v0d.im, v0d.im, temp1.im);
    mpz_sub(v0d.re, v0d.re, this->v0.re);
    mpz_sub(v0d.im, v0d.im, this->v0.im);

    ExtendedMumford ret(f, h, u1d, u0d, v1d, v0d, U1d, U0d);
    return ret;
}


ExtendedMumford ExtendedMumford::zero(const ExtendedPolynomial& f, const ExtendedPolynomial& h){
    ExtendedMumford zero(f, h);
    return zero;
}

bool ExtendedMumford::isZero() const{
    if(v1.isZero() && v0.isZero()){
        if(u1.isZero() && u0 == ExtendedNumber::ONE()){
            return true;
        }
    }
    return false;
}

void ExtendedMumford::print() const{
    std::cerr << "[" << u1 << ", " << u0 << "]" << std::endl;
    std::cerr << "[" << v1 << ", " << v0 << "]" << std::endl;
    return;
}


