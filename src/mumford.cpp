#include "mumford.hpp"

Mumford::Mumford(){
    this->f = Polynomial();
    this->h = Polynomial();
    this->u2 = Number::ONE();
    this->u1 = Number::ZERO();
    this->u0 = Number::ONE();
    this->v1 = Number::ZERO();
    this->v0 = Number::ZERO();
    this->U1 = Number::ZERO();
    this->U0 = Number::ZERO();
}

Mumford::Mumford(Polynomial f, Polynomial h){
    this->f = f;
    this->h = h;
    this->u2 = Number::ONE();
    this->u1 = Number::ZERO();
    this->u0 = Number::ONE();
    this->v1 = Number::ZERO();
    this->v0 = Number::ZERO();
    this->U1 = Number::ZERO();
    this->U0 = Number::ZERO();    
}

Mumford::Mumford(Polynomial f, Polynomial h, Number u1, Number u0, Number v1, Number v0){
    this->f = f;
    this->h = h;
    this->u2 = Number::ONE();
    this->u1 = u1;
    this->u0 = u0;
    this->v1 = v1;
    this->v0 = v0;
    this->U1 = u1 * u1;
    this->U0 = u1 * u0;
}

Mumford::Mumford(Polynomial f, Polynomial h, Number u1, Number u0, Number v1, Number v0, Number U1, Number U0){
    this->f = f;
    this->h = h;
    this->u2 = Number::ONE();
    this->u1 = u1;
    this->u0 = u0;
    this->v1 = v1;
    this->v0 = v0;
    this->U1 = U1;
    this->U0 = U0;
}

// 被約因子 D を受け取り，対応する Mumford 表現を返す．
Mumford::Mumford(Polynomial f, Polynomial h, Divisor d){
    this->f = f;
    this->h = h;

    int count = d.points.size();
    if(count == 0){
        this->u2 = Number::ZERO();
        this->u1 = Number::ZERO();
        this->u0 = Number::ONE();
        this->v1 = Number::ZERO();
        this->v0 = Number::ZERO();
    }else if(count == 1){
        auto p = d.points[0];
        int multiplicity = p.second;
        if(multiplicity == 1){
            this->u1 = Number::ONE();
            this->u0 = -p.first.x;
            this->v1 = Number::ZERO();
            this->v0 = p.first.y;
        }else{
            // TODO : 重複度が 2 の場合．
        }
    }else if(count == 2){
        // 2 つの異なる点が含まれる場合．
        Point p1 = d.points[0].first;
        Point p2 = d.points[1].first;

        this->u2 = Number::ONE();
        this->u1 = (p1.x * p2.x);
        this->u0 = -(p1.x + p2.x);

        Number c1 = (p1.y - p2.y) / (p1.x - p2.x);
        Number c0 = (p1.x * p2.y - p2.x * p1.y) / (p1.x - p2.x);
        this->v1 = c1;
        this->v0 = c0;
    }
}

Mumford Mumford::CostelloScalarMultiple (const mpz_class& k_) const{
    Polynomial f = this->f;
    Polynomial h = this->h;
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

    Mumford D = Mumford::zero(f, h);
    Mumford now = *this;
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

Mumford Mumford::operator * (const mpz_class& k_) const{
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

    Mumford D = Mumford::zero(f, h);
    Mumford now = *this;
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

Mumford Mumford::operator + (const Mumford& m) const{
    if(this->f.deg == 6){
        return this->LangeAdd(m);
    }else{
        return this->CantorAdd(m);
    }
}

Mumford Mumford::CostelloAdd(const Mumford& m) const{
    //std::cerr << "Costello Addition." << std::endl;
    Number u11 = this->u1;
    Number u10 = this->u0;
    Number u21 = m.u1;
    Number u20 = m.u0;

    Number v11 = this->v1;
    Number v10 = this->v0;
    Number v21 = m.v1;
    Number v20 = m.v0;

    Number f6;
    if(this->f.deg == 6){
        f6 = this->f.coeff[6];
    }else{
        f6 = Number::ZERO();
    }
    Number f5 = this->f.coeff[5];
    Number f4 = this->f.coeff[4];

    Number U11 = this->U1;
    Number U21 = m.U1;
    Number U10 = this->U0;
    Number U20 = m.U0;

    Number temp1;

    Number u1S, v0D, v1D;
    mpz_add(u1S.value, u11.value, u21.value);
    mpz_sub(v0D.value, v10.value, v20.value);
    mpz_sub(v1D.value, v11.value, v21.value);

    Number M1, M2, M3, M4;
    mpz_sub(M1.value, U11.value, U21.value);
    mpz_sub(M1.value, M1.value, u10.value);
    mpz_add(M1.value, M1.value, u20.value);
    mpz_sub(M2.value, U20.value, U10.value);
    mpz_sub(M3.value, u11.value, u21.value);
    mpz_sub(M4.value, u20.value, u10.value);

    Number t1, t2, t3, t4;
    mpz_sub(t1.value, M2.value, v0D.value);
    mpz_sub(temp1.value, v1D.value, M1.value);
    mpz_mul(t1.value, t1.value, temp1.value);
    mpz_mod(t1.value, t1.value, Number::CHARA);

    mpz_add(t2.value, M2.value, v0D.value);
    mpz_add(temp1.value, v1D.value, M1.value);
    mpz_mul(t2.value, t2.value, temp1.value);
    mpz_mod(t2.value, t2.value, Number::CHARA);
    mpz_neg(t2.value, t2.value);

    mpz_sub(t3.value, M4.value, v0D.value);
    mpz_sub(temp1.value, v1D.value, M3.value);
    mpz_mul(t3.value, t3.value, temp1.value);
    mpz_mod(t3.value, t3.value, Number::CHARA);

    mpz_add(t4.value, M4.value, v0D.value);
    mpz_add(temp1.value, v1D.value, M3.value);
    mpz_mul(t4.value, t4.value, temp1.value);
    mpz_mod(t4.value, t4.value, Number::CHARA);
    mpz_neg(t4.value, t4.value);

    Number l2_num, l3_num;
    mpz_sub(l2_num.value, t1.value, t2.value);
    mpz_sub(l3_num.value, t3.value, t4.value);

    Number d;
    mpz_sub(d.value, M4.value, M2.value);
    mpz_add(temp1.value, M1.value, M3.value);
    mpz_mul(d.value, d.value, temp1.value);
    mpz_mod(d.value, d.value, Number::CHARA);
    mpz_add(d.value, d.value, d.value);
    mpz_add(d.value, d.value, t3.value);
    mpz_add(d.value, d.value, t4.value);
    mpz_sub(d.value, d.value, t1.value);
    mpz_sub(d.value, d.value, t2.value);

    Number A, B, C, d_inv, d_shifted_inv;
    mpz_mul(A.value, d.value, d.value);
    mpz_mod(A.value, A.value, Number::CHARA);
    mpz_mul(temp1.value, f6.value, A.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);    
    mpz_mul(B.value, l3_num.value, l3_num.value);
    mpz_mod(B.value, B.value, Number::CHARA);
    mpz_sub(B.value, B.value, temp1.value);
    mpz_mul(C.value, d.value, B.value);
    mpz_mod(C.value, C.value, Number::CHARA);
    mpz_invert(C.value, C.value, Number::CHARA);
    mpz_mul(d_inv.value, B.value, C.value);
    mpz_mod(d_inv.value, d_inv.value, Number::CHARA);
    mpz_mul(d_shifted_inv.value, d.value, A.value);
    mpz_mod(d_shifted_inv.value, d_shifted_inv.value, Number::CHARA);
    mpz_mul(d_shifted_inv.value, d_shifted_inv.value, C.value);
    mpz_mod(d_shifted_inv.value, d_shifted_inv.value, Number::CHARA);

    Number l2, l3;
    mpz_mul(l2.value, l2_num.value, d_inv.value);
    mpz_mod(l2.value, l2.value, Number::CHARA);
    mpz_mul(l3.value, l3_num.value, d_inv.value);
    mpz_mod(l3.value, l3.value, Number::CHARA);

    mpz_mul(temp1.value, l2.value, l3.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);

    Number u1dd;
    mpz_sub(u1dd.value, f5.value, temp1.value);
    mpz_sub(u1dd.value, u1dd.value, temp1.value);
    mpz_mul(u1dd.value, u1dd.value, d_shifted_inv.value);
    mpz_mod(u1dd.value, u1dd.value, Number::CHARA);
    mpz_add(u1dd.value, u1dd.value, u1S.value);
    mpz_neg(u1dd.value, u1dd.value);

    Number u0dd;
    mpz_sub(u0dd.value, u10.value, U11.value);
    mpz_mul(u0dd.value, u0dd.value, l3.value);
    mpz_mod(u0dd.value, u0dd.value, Number::CHARA);
    mpz_mul(temp1.value, l2.value, u11.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);    
    mpz_add(u0dd.value, u0dd.value, temp1.value);
    mpz_add(u0dd.value, u0dd.value, v11.value);
    mpz_mul(u0dd.value, u0dd.value, l3.value);
    mpz_mod(u0dd.value, u0dd.value, Number::CHARA);
    mpz_add(u0dd.value, u0dd.value, u0dd.value);
    mpz_mul(temp1.value, l2.value, l2.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);    
    mpz_add(u0dd.value, u0dd.value, temp1.value);
    mpz_sub(u0dd.value, u0dd.value, f4.value);
    mpz_mul(u0dd.value, u0dd.value, d_shifted_inv.value);
    mpz_mod(u0dd.value, u0dd.value, Number::CHARA);
    mpz_mul(temp1.value, u21.value, u11.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);    
    mpz_sub(u0dd.value, u0dd.value, temp1.value);
    mpz_sub(u0dd.value, u0dd.value, u10.value);
    mpz_sub(u0dd.value, u0dd.value, u20.value);
    mpz_mul(temp1.value, u1S.value, u1dd.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);    
    mpz_sub(u0dd.value, u0dd.value, temp1.value);

    Number U1dd, U0dd;
    mpz_mul(U1dd.value, u1dd.value, u1dd.value);
    mpz_mod(U1dd.value, U1dd.value, Number::CHARA);
    mpz_mul(U0dd.value, u1dd.value, u0dd.value);
    mpz_mod(U0dd.value, U0dd.value, Number::CHARA);

    Number v1dd, v0dd;
    mpz_sub(temp1.value, u1dd.value, u11.value);
    mpz_mul(temp1.value, l2.value, temp1.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);

    mpz_sub(v1dd.value, u0dd.value, U1dd.value);
    mpz_add(v1dd.value, v1dd.value, U11.value);
    mpz_sub(v1dd.value, v1dd.value, u10.value);
    mpz_mul(v1dd.value, v1dd.value, l3.value);
    mpz_mod(v1dd.value, v1dd.value, Number::CHARA);
    mpz_add(v1dd.value, v1dd.value, temp1.value);
    mpz_sub(v1dd.value, v1dd.value, v11.value);

    mpz_sub(temp1.value, u0dd.value, u10.value);
    mpz_mul(temp1.value, l2.value, temp1.value);
    mpz_mod(temp1.value, temp1.value, Number::CHARA);

    mpz_sub(v0dd.value, U10.value, U0dd.value);
    mpz_mul(v0dd.value, v0dd.value, l3.value);
    mpz_mod(v0dd.value, v0dd.value, Number::CHARA);
    mpz_add(v0dd.value, v0dd.value, temp1.value);
    mpz_sub(v0dd.value, v0dd.value, v10.value);

    return Mumford(f, h, u1dd, u0dd, v1dd, v0dd, U1dd, U0dd);
}

Mumford Mumford::HarleyAdd(const Mumford& m) const{
    std::cerr << "Harley Addition." << std::endl;
    Number u11 = this->u1;
    Number u10 = this->u0;
    Number u21 = m.u1;
    Number u20 = m.u0;

    Number v11 = this->v1;
    Number v10 = this->v0;
    Number v21 = m.v1;
    Number v20 = m.v0;

    Number h0 = this->h.coeff[0];
    Number h1 = this->h.coeff[1];
    Number h2 = this->h.coeff[2];

    Number f4 = this->f.coeff[4];

    // 1. u1, u2 の終結式を計算．
    Number z1 = u11 - u21;
    Number z2 = u20 - u10;
    Number z3 = u11 * z1 + z2;
    Number r = z2 * z3 + z1 * z1 * u10;

    // 2. u2 の almost inverse (mod u1) を計算．
    Number inv1 = z1;
    Number inv0 = z3;

    // 3. s' を計算．
    Number w0 = v10 - v20;
    Number w1 = v11 - v21;
    Number w2 = inv0 * w0;
    Number w3 = inv1 * w1;
    Number s1d = (inv0 + inv1) * (w0 + w1) - w2 - w3 * (Number::ONE() + u11);
    Number s0d = w2 - u10 * w3;

    if(!s1d.isZero()){
        // 4. s'' を計算．
        Number w1 = (r * s1d).inv();
        Number w2 = r * w1;
        Number w3 = s1d * s1d * w1;
        Number w4 = r * w2;
        Number w5 = w4 * w4;
        Number s0dd = s0d * w2;

        // 5. l' を計算．
        Number l2d = u21 + s0dd;
        Number l1d = u21 * s0dd + u20;
        Number l0d = u20 * s0dd;


        // 6. u' を計算．
        Number u0d = (s0dd - u11) * (s0dd - z1 + h2 * w4) - u10 + l1d + (h1 + v21 * 2) * w4 + (u21 * 2 + z1 - f4) * w5;
        Number u1d = s0dd * 2 - z1 + h2 * w4 - w5;

        // 7. v' を計算．
        w1 = l2d - u1d;
        w2 = u1d * w1 + u0d - l1d;
        Number v1d = w2 * w3 - v21 - h1 + h2 * u1d;
        w2 = u0d * w1 - l0d;
        Number v0d = w2 * w3 - v20 - h0 + h2 * u0d;

        Mumford ret(f, h, u1d, u0d, v1d, v0d);
        return ret;
    }else{
        //std::cerr << "Special case." << std::endl;
        // サブルーチン

        // 4'. s を計算．
        Number inv = r.inv();
        Number s0 = s0d * inv;

        // 5'. u' を計算．
        Number u0d = f4 - u21 - u11 - s0*s0 - s0 * h2;

        // 6'. v' を計算．
        Number w1 = s0 * (u21 + u0d) + h1 + v21 + h2 * u0d;
        Number w2 = s0 + v20 + h0;
        Number v0d = u0d * w1 - w2;

        Mumford ret(f, h, Number::ONE(), u0d, Number::ZERO(), v0d);
        return ret;
    }
}

Mumford Mumford::HarleyAddDegenerated(const Mumford& m) const{
    Number u11 = this->u1;
    Number u10 = this->u0;
    Number u21 = m.u1;
    Number u20 = m.u0;

    Number v11 = this->v1;
    Number v10 = this->v0;
    Number v21 = m.v1;
    Number v20 = m.v0;

    Number h0 = this->h.coeff[0];
    Number h1 = this->h.coeff[1];
    Number h2 = this->h.coeff[2];

    Number f3 = this->f.coeff[3];
    Number f4 = this->f.coeff[4];

    // 1. r を計算．
    Number r = u20 - (u21 - u10) * u10;

    // 2. u2 の almost inverse (mod u1) を計算．
    Number inv = r.inv();

    // 3. s を計算．
    Number s0 = inv * (v10 - v20 + v21 * u10);

    // 4. l を計算．
    Number l1 = s0 * u21;
    Number l0 = s0 * u20;

    // 5. k を計算．
    Number k2 = f4 - u21;
    Number k1 = f3 - (f4 - u21) * u21 - v21 * h2 - u20;

    // 6. u' を計算．
    Number u1d = k2 - s0 * s0 - s0 * h2 - u10;
    Number u0d = k1 - s0 * (l1 + h1 + v21 * 2) - u10 * u1d;

    // 7. v' を計算．
    Number v1d = (h2 + s0) * u1d - (h1 + l1 + v21);
    Number v0d = (h2 + s0) * u0d - (h0 + l0 + v20);

    Mumford ret(f, h, u1d, u0d, v1d, v0d);
    return ret;
}

Mumford Mumford::LangeAdd(const Mumford& m) const{
    //std::cerr << "Lange Addition." << std::endl;
    Number u11 = this->u1;
    Number u10 = this->u0;
    Number u21 = m.u1;
    Number u20 = m.u0;

    Number v11 = this->v1;
    Number v10 = this->v0;
    Number v21 = m.v1;
    Number v20 = m.v0;

    Number f6 = this->f.coeff[6];
    Number f5 = this->f.coeff[5];
    Number f4 = this->f.coeff[4];

    // 1. u1, u2 の終結式を計算．
    Number z1 = u11 - u21;
    Number z2 = u20 - u10;
    Number z3 = u11 * z1 + z2;
    Number r = z2 * z3 + z1 * z1 * u10;

    // 2. u2 の almost inverse (mod u1) を計算．
    Number inv1 = z1;
    Number inv0 = z3;

    // 3. s' を計算．
    Number w0 = v10 - v20;
    Number w1 = v11 - v21;
    Number w2 = inv0 * w0;
    Number w3 = inv1 * w1;
    Number s1d = (inv0 + inv1) * (w0 + w1) - w2 - w3 * (Number::ONE() + u11);
    Number s0d = w2 - u10 * w3;

    if(s1d.isZero()){
        //std::cerr << "Special case." << std::endl;
        // todo: ここの場合分けを厳密に書く．
        Mumford ret(f, h);
        return ret;
    }

    // 4. l' を計算．
    Number l3d = s1d;
    Number l2d = u21 * s1d;
    Number l0d = u20 * s0d;
    Number l1d = (s1d + s0d) * (u21 + u20) - l2d - l0d; //u21 * s0d + u20 * s1d; 
    l2d = l2d + s0d;

    // 5. u' を計算．
    Number k4 = f6;
    Number k3 = f5 - f6 * u21;
    Number k2 = f4 - f6 * u20 - (f5 - f6 * u21) * u21;

    Number t4 = s1d * l3d - k4 * r * r;
    Number t3 = s1d * l2d + s0d * l3d - k3 * r * r;
    Number t2 = r * v21;
    t2 = t2 + t2;
    t2 = s1d * (t2 + l1d) + s0d * l2d - k2 * r * r;

    Number u0d = t2 - t4 * u10 - (t3 - t4 * u11) * u11;
    Number u1d = t3 - t4 * u11;
    Number u2d = t4;

    // 6. u' をモニックにする．
    w1 = (u2d * r).inv();
    w2 = w1 * r;
    w3 = w1 * u2d;
    u1d = u1d * w2;
    u0d = u0d * w2;
    // ud の計算まで正しい．

    // 7. v' を計算．
    Number v1d = (-l1d + (u0d - u1d * u1d) * l3d + u1d * l2d) * w3 - v21;
    Number v0d = (-l0d  - u1d * u0d * l3d + u0d * l2d) * w3 - v20;

    Mumford ret(f, h, u1d, u0d, v1d, v0d);

    return ret;
}

Mumford Mumford::LangeDoubling() const{
    //std::cerr << "Lange Doubling." << std::endl;
    Number u1 = this->u1;
    Number u0 = this->u0;
    Number v1 = this->v1;
    Number v0 = this->v0;

    Number f2 = this->f.coeff[2];
    Number f3 = this->f.coeff[3];
    Number f4 = this->f.coeff[4];
    Number f5 = this->f.coeff[5];
    Number f6 = this->f.coeff[6];

    // 44M, 6S, I

    // 1. v~ (=2v) と u の終結式を計算．
    // 3M, 2S
    Number v1t = v1 + v1;
    Number v0t = v0 + v0;

    Number w0 = v1 * v1;
    Number u1s = u1 * u1;
    Number w1 = u1s;
    Number w2 = w0 + w0;
    w2 = w2 + w2;
    Number w3 = u1 * v1t;
    Number r = u0 * w2 + v0t * (v0t - w3);

    // 2. r の almost inverse を計算．
    Number inv1d = -v1t;
    Number inv0d = v0t - w3;

    // 3. k を計算．
    // 11M
    Number k4 = f6;
    Number k4u0 = k4 * u0;
    Number k4u1 = k4 * u1;
    Number k3 = f5 - k4u1;
    Number k3u0 = k3 * u0;
    Number k2 = f4 - k4u0 - k3 * u1;
    Number k1 = f3 - k3u0 - k2 * u1;
    Number k0 = f2 - w0 - k2 * u0 - k1 * u1;
    Number u1kd = u1 * (k3 - k4u1);
    Number k1d = k1 + w1 * (k3 - k4u1) - k3u0 + u1 * (k4u0 + k4u0 - k2);
    Number k0d = k0 + u0 * (u1kd + k4u0 - k2);

    // 4. s' を計算．
    // 5M
    w0 = k0d * inv0d;
    w1 = k1d * inv1d;
    Number s1d = (inv0d + inv1d) * (k0d + k1d) - w0 - w1 * (Number::ONE() + u1);
    Number s0d = w0 - u0 * w1;

    if(s1d.isZero()){
        std::cerr << "Special case." << std::endl;
        // サブルーチン
        Mumford ret(f, h);
        return ret;
    }

    // 5. u' を計算．
    // 8M, 3S
    Number rs = r * r;
    Number u2d = s1d * s1d - rs * f6;
    Number u1d = s1d * s0d;
    u1d = u1d + u1d;
    u1d = u1d - rs * (f5 - k4u1 - k4u1);
    Number v1s1dr_2 = v1 * s1d * r;
    v1s1dr_2 = v1s1dr_2 + v1s1dr_2;
    Number u0d = v1 * s1d * r;
    u0d = u0d + u0d;
    u0d = u0d + s0d * s0d;
    Number ft = (f5 - k4u1 - k4u1);
    ft = ft + ft;
    u0d = u0d - rs * (f4 - (u0 + u0 + u1s) * f6 - u1 * ft);

    // 6. r と u2d の逆元を計算．
    // 3M, I
    w0 = (r * u2d).inv();
    w1 = w0 * r;
    w2 = w0 * u2d;

    // 7. u' を計算．
    // 2M
    u1d = u1d * w1;
    u0d = u0d * w1;

    // 8. v' を計算．
    // 12M, 1S
    Number l3 = s1d * w2;
    Number l2 = (s1d * u1 + s0d) * w2;
    Number l1 = (s1d * u0 + s0d * u1) * w2;
    Number l0 = s0d * u0 * w2;
    Number v1d = l3 * (u1d * u1d - u0d) - l2 * u1d + l1 + v1;
    Number v0d = (l3 * u1d - l2) * u0d + l0 + v0;

    Mumford ret(f, h, u1d, u0d, -v1d, -v0d);
    return ret;
}

Mumford Mumford::CostelloDoubling() const{
    //std::cerr << "Costello Doubling." << std::endl;
    Number u1 = this->u1;
    Number u0 = this->u0;
    Number v1 = this->v1;
    Number v0 = this->v0;

    Number U1 = this->U1;
    Number U0 = this->U0;

    Number f2 = this->f.coeff[2];
    Number f3 = this->f.coeff[3];
    Number f4 = this->f.coeff[4];
    Number f5 = this->f.coeff[5];
    Number f6 = this->f.coeff[6];

    // 32M, 6S, I

    Number vv = v1 * v1;
    Number va = (v1 + u1) * (v1 + u1) - vv - U1;

    Number M1 = (v0 - va);
    M1 = M1 + M1;
    Number M2 = (U1 + U1 + u0) * v1;
    M2 = M2 + M2;
    Number M3 = -(v1 + v1);
    Number M4 = va + v0 + v0;

    Number f4_2 = f4 + f4;
    Number f6u0 = f6 * u0;
    Number f6U1 = f6 * U1;
    Number f5u0 = f5 * u0;
    Number f5u1 = f5 * u1;
    Number f5u1_2 = f5u1 + f5u1;
    Number f6u0_3 = f6u0 + f6u0 + f6u0;
    Number f6u0_6 = f6u0_3 + f6u0_3;
    Number f6U1_2 = f6U1 + f6U1;
    Number f6U1_4 = f6U1_2 + f6U1_2;

    Number z11 = f5u1_2 - f6U1_2 - f6U1 - f4;
    z11 = z11 * U1;
    Number z12 = f5u1_2 + f6u0_3 - f4_2;
    z12 = z12 * u0;

    Number z1 = z11 + z12 - vv + f2;
    Number z2 = (f6u0_6 - f6U1_4 + f5u1_2 + f5u1 - f4_2) * u1 - f5u0 - f5u0 + f3;

    Number t11 = M2 - z1;
    Number t12 = z2 - M1;
    Number t21 = -(z1 + M2);
    Number t22 = z2 + M1;
    Number t31 = -z1 + M4;
    Number t32 = z2 - M3;
    Number t41 = -(z1 + M4);
    Number t42 = z2 + M3;

    Number t1 = t11 * t12;
    Number t2 = t21 * t22;
    Number t3 = t31 * t32;
    Number t4 = t41 * t42;

    Number l2_num = t1 - t2;
    Number l3_num = t3 - t4;

    Number d11 = M4 - M2;
    Number d12 = M1 + M3;
    Number d = d11 * d12;
    d = d + d;
    d = d + t3 + t4 - t1 - t2;

    Number A = d * d;
    Number B1 = l3_num * l3_num;
    Number B2 = f6 * A;
    Number B = B1 - B2;
    Number C = (d * B).inv();

    Number d_inv = B * C;
    Number d_shifted_inv = A * C * d;

    Number l2 = l2_num * d_inv;
    Number l3 = l3_num * d_inv;

    Number l2l3 = l2 * l3;
    Number u1d = (l2l3 + l2l3 - f5) * d_shifted_inv - u1 - u1;
    Number u0d = ((u0 - U1) * l3 + l2 * u1 + v1) * l3;
    u0d = u0d + u0d;
    u0d = u0d + l2 * l2 - f4;
    u0d = u0d * d_shifted_inv;
    Number u1u1d_2 = u1 * u1d;
    u1u1d_2 = u1u1d_2 + u1u1d_2;
    u0d = u0d - u0 - u0 - u1u1d_2 - U1;
    Number U1d = u1d * u1d;
    Number U0d = u1d * u0d;

    Number v1d = (u0d - U1d + U1 - u0) * l3 + (u1d - u1) * l2 - v1;
    Number v0d = (U0 - U0d) * l3 + (u0d - u0) * l2 - v0;

    Mumford ret(f, h, u1d, u0d, v1d, v0d, U1d, U0d);
    return ret;
}

Mumford Mumford::CantorAdd(const Mumford& m) const{
    std::cerr << "Cantor Addition." << std::endl;

    Polynomial h = this->h;

    Polynomial u1 = Polynomial(2);
    u1.coeff[2] = Number::ONE();
    u1.coeff[1] = this->u1;
    u1.coeff[0] = this->u0;
    Polynomial v1 = Polynomial(1, this->v0, this->v1);
    Polynomial u2 = Polynomial(2);
    u2.coeff[2] = Number::ONE();
    u2.coeff[1] = m.u1;
    u2.coeff[0] = m.u0;
    Polynomial v2 = Polynomial(1, m.v0, m.v1);

    auto tup1 = Polynomial::extended_gcd(u1, u2);
    Polynomial d1 = std::get<0>(tup1);
    Polynomial e1 = std::get<1>(tup1);
    Polynomial e2 = std::get<2>(tup1);

    auto tup2 = Polynomial::extended_gcd(d1, v1 + v2 + h);
    Polynomial d = std::get<0>(tup2);
    Polynomial c1 = std::get<1>(tup2);
    Polynomial c2 = std::get<2>(tup2);

    Polynomial s1 = c1 * e1;
    Polynomial s2 = c1 * e2;

    Polynomial u = u1 * u2 / (d * d);
    Polynomial v = ((s1 * u1 * v2 + s2 * u2 * v1 + c2 * (v1 * v2 + f)) / d) % u;

    while(u.deg > GENUS){
        Polynomial ud = (f - (v * h) - (v * v)) / u;
        Polynomial vr = (-(h + v)) % ud;

        u = ud;
        v = vr;
    }

    Number lc = u.coeff[u.deg];
    u = u * lc.inv();

    Number v1n = (v.deg == 1) ? v.coeff[1] : Number::ZERO();

    Mumford ret(f, h, u.coeff[1], u.coeff[0], v1n, v.coeff[0]);

    return ret;
}

Mumford Mumford::inv(){
    Polynomial f = this->f;
    Polynomial h = this->h;
    Polynomial u(2);
    u.coeff[0] = u0;
    u.coeff[1] = u1;
    u.coeff[2] = u2;
    Polynomial v(1, v0, v1);

    Polynomial vd = (-(h + v)) % u;

    Mumford inv(f, h, u1, u0, vd.coeff[1], vd.coeff[0]);
    return inv;
}

Mumford Mumford::zero() const{
    Polynomial f = this->f;
    Polynomial h = this->h;

    Mumford zero(f, h);
    return zero;
}

Mumford Mumford::zero(const Polynomial& f, const Polynomial& h){
    Mumford zero(f, h);
    return zero;
}

void Mumford::print() const{
    Polynomial u(2);
    u.coeff[0] = u0;
    u.coeff[1] = u1;
    u.coeff[2] = Number::ONE();
    Polynomial v(1, v1, v0);
    std::cerr << "[" << u << ", " << v << "]" << std::endl;
    return;
}

bool Mumford::isZero() const{
    if(v1.isZero() && v0.isZero()){
        if(u1.isZero() && u0 == Number::ONE()){
            return true;
        }
    }
    return false;
}
