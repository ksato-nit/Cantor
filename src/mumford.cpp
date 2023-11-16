#include "mumford.hpp"

Mumford::Mumford(){
    this->f = Polynomial();
    this->h = Polynomial();
    this->u = Polynomial();
    this->v = Polynomial();
}

Mumford::Mumford(Polynomial f, Polynomial h){
    this->f = f;
    this->h = h;
    this->u = Polynomial();
    this->v = Polynomial();
}

Mumford::Mumford(Polynomial f, Polynomial h, Polynomial u, Polynomial v){
    this->f = f;
    this->h = h;
    this->u = u;
    this->v = v;
}

// 被約因子 D を受け取り，対応する Mumford 表現を返す．
Mumford::Mumford(Polynomial f, Polynomial h, Divisor d){
    this->f = f;
    this->h = h;

    int count = d.points.size();
    if(count == 0){
        Polynomial ONE(0, 1);
        Polynomial ZERO(0, 0);
        this->u = ONE;
        this->v = ZERO;
    }else if(count == 1){
        auto p = d.points[0];
        int multiplicity = p.second;
        if(multiplicity == 1){
            Polynomial u(1, Number::ONE(), -p.first.x);
            Polynomial v(0, p.first.y);
            this->u = u;
            this->v = v;
        }else{
            // TODO : 重複度が 2 の場合．
        }
    }else if(count == 2){
        // 2 つの異なる点が含まれる場合．
        Point p1 = d.points[0].first;
        Point p2 = d.points[1].first;

        std::vector<Number> uc = {(p1.x * p2.x), -(p1.x + p2.x), Number::ONE()};
        Polynomial u(uc);
        this->u = u;

        Number c1 = (p1.y - p2.y) / (p1.x - p2.x);
        Number c0 = (p1.x * p2.y - p2.x * p1.y) / (p1.x - p2.x);
        Polynomial v(1, c0, c1);
        this->v = v;

    }
}

Mumford Mumford::operator + (const Mumford& m) const{
    Polynomial u1 = this->u;
    Polynomial u2 = m.u;

    Polynomial v1 = this->v;
    Polynomial v2 = m.v;

    Number u10 = u1.coeff[0];
    Number u20 = u2.coeff[0];
    Number v10 = v1.coeff[0];
    Number v20 = v2.coeff[0];

    if(u1.deg > u2.deg){
        //std::cerr << "Flip." << std::endl;
        return m + *this;
    }

    if(u1.deg == 0){
        return m;
    }

    /*
    if(u1.deg == 1){
        if(u2.deg == 1){
            if(u1 == u2){
                Number h_eval = h.eval(-u10);
                Polynomial h_eval_as_poly(0, h_eval);
                if(-v1\ == v2 + h_eval_as_poly){
                    return Mumford::zero();
                }else{
                    Polynomial u = u1 * u1;
                    Polynomial v; // TODO : Lange p. 6 をもとに書く．
                    Mumford ret(f, h, u, v);
                }
            }else{
                Polynomial u = u1 * u2;

                Number c1 = (v20 - v10) / (u10 - u20);
                Number c0 = (v20 * u10 - v10 * u20) / (u10 - u20);
                Polynomial v(1, c0, c1);
            }
        }else{
            Number u10 = u1.coeff[0];
            Number u21 = u2.coeff[1];
            Number h_eval = h.eval(-u10);
            Number v2_eval = v2.eval(-u10);

            if(v2_eval == v10 + h_eval){
                Polynomial u(1, Number::ONE(), u21 - u10);
                Polynomial v(0, v20 * (u10 - u21));
                Mumford ret(f, h, u, v);
                return ret;
            }else{
                // 2. (b) ii 後半
            }
            //std::cerr << "Degenerated." << std::endl;
            Mumford ret = this->HarleyAddDegenerated(m);
            return ret;
        }
    }

    if(u1 == u2 && v1 == v2){
        Mumford ret = this->doubling();
        return ret;
    }else{
        Mumford ret = this->HarleyAdd(m);
        return ret;
    }
    */

    // deg u1 = deg u2 = 2
    if(this->f.deg == 6){
        return this->LangeAdd(m);
    }else{
        return this->CantorAdd(m);
    }
}

Mumford Mumford::CostelloAdd(const Mumford& m) const{
    std::cout << "Costello Addition." << std::endl;
    Polynomial u1 = this->u;
    Polynomial v1 = this->v;
    Polynomial u2 = m.u;
    Polynomial v2 = m.v;

    Number u11 = u1.coeff[1];
    Number u10 = u1.coeff[0];
    Number u21 = u2.coeff[1];
    Number u20 = u2.coeff[0];

    Number v11 = v1.coeff[1];
    Number v10 = v1.coeff[0];
    Number v21 = v2.coeff[1];
    Number v20 = v2.coeff[0];

    Number f6;
    if(this->f.deg == 6){
        f6 = this->f.coeff[6];
    }else{
        f6 = Number::ZERO();
    }
    Number f5 = this->f.coeff[5];
    Number f4 = this->f.coeff[4];

    // TODO: U1, U0 を座標に含めて保持するようにする．
    Number U11 = u11 * u11;
    Number U21 = u21 * u21;
    Number U10 = u11 * u10;
    Number U20 = u21 * u20;

    Number u1S = u11 + u21;
    Number v0D = v10 - v20;
    Number v1D = v11 - v21;

    Number M1 = U11 - U21 - u10 + u20;
    Number M2 = U20 - U10;
    Number M3 = u11 - u21;
    Number M4 = u20 - u10;

    //std::cout << M1 << " " << M2 << " " << M3 << " " << M4 << std::endl;

    Number t1 = (M2 - v0D) * (v1D - M1);
    Number t2 = -(v0D + M2) * (v1D + M1);
    Number t3 = (M4 - v0D) * (v1D - M3);
    Number t4 = -(v0D + M4) * (v1D + M3);

    //std::cout << t1 << " " << t2 << " " << t3 << " " << t4 << std::endl;

    Number l2_num = t1 - t2;
    Number l3_num = t3 - t4;

    Number d = (M4 - M2) * (M1 + M3);
    d = d + d;
    d = d + t3 + t4 - t1 - t2;

    Number A = d * d;
    Number B = l3_num * l3_num - f6 * A;
    Number C = (d * B).inv();
    Number d_inv = B * C;
    Number d_shifted_inv = d * A * C;

    Number l2 = l2_num * d_inv;
    Number l3 = l3_num * d_inv;

    /*
    std::cout << "d : " << d << std::endl;
    std::cout << "B : " << B << std::endl;
    std::cout << "d_shifted_inv : " << d_shifted_inv << std::endl;
    
    */
    std::cout << l2_num << " " << l3_num << " " << d << std::endl;
    //std::cout << l3 << " " << l2 << " " << std::endl;

    Number l0 = v10 + l2 * u10 - l3 * U10;
    Number l1 = v11 + l2 * u11 - l3 * (U11 - u10);

    Number u1dd = -(u1S + (f5 - l2 * l3 - l2 * l3) * d_shifted_inv);

    Number u0dd = l3 * (l3 * (u10 - U11) + l2 * u11 + v11);
    u0dd = u0dd + u0dd;
    u0dd = u0dd + l2 * l2 - f4;
    u0dd = u0dd * d_shifted_inv;
    u0dd = u0dd - u11 * u21 - u10 - u20 - u1S * u1dd;

    Number U1dd = u1dd * u1dd;
    Number U0dd = u1dd * u0dd;

    Number v1dd = l3 * (u0dd - U1dd + U11 - u10) + l2 * (u1dd - u11) - v11;
    Number v0dd = l3 * (U10 - U0dd) + l2 * (u0dd - u10) - v10;

    Polynomial u(2);
    u.coeff[2] = Number::ONE();
    u.coeff[1] = u1dd;
    u.coeff[0] = u0dd;

    Polynomial v(1, v0dd, v1dd);

    return Mumford(f, h, u, v);
}

Mumford Mumford::HarleyAdd(const Mumford& m) const{
    std::cout << "Harley Addition." << std::endl;
    Polynomial u1 = this->u;
    Polynomial v1 = this->v;
    Polynomial u2 = m.u;
    Polynomial v2 = m.v;

    Number u11 = u1.coeff[1];
    Number u10 = u1.coeff[0];
    Number u21 = u2.coeff[1];
    Number u20 = u2.coeff[0];

    Number v11 = v1.coeff[1];
    Number v10 = v1.coeff[0];
    Number v21 = v2.coeff[1];
    Number v20 = v2.coeff[0];

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

        Polynomial u(2);
        Polynomial v(1);

        u.coeff[2] = Number::ONE();
        u.coeff[1] = u1d;
        u.coeff[0] = u0d;

        v.coeff[1] = v1d;
        v.coeff[0] = v0d;

        Mumford ret(f, h, u, v);
        return ret;
    }else{
        std::cout << "Special case." << std::endl;
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

        Polynomial u(1);
        Polynomial v(0);

        u.coeff[1] = Number::ONE();
        u.coeff[0] = u0d;

        v.coeff[0] = v0d;

        Mumford ret(f, h, u, v);
        return ret;
    }
}

Mumford Mumford::HarleyAddDegenerated(const Mumford& m) const{
    Polynomial u1 = this->u;
    Polynomial v1 = this->v;
    Polynomial u2 = m.u;
    Polynomial v2 = m.v;

    Number u11 = u1.coeff[1];
    Number u10 = u1.coeff[0];
    Number u21 = u2.coeff[1];
    Number u20 = u2.coeff[0];

    Number v11 = v1.coeff[1];
    Number v10 = v1.coeff[0];
    Number v21 = v2.coeff[1];
    Number v20 = v2.coeff[0];

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

    Polynomial u(2);
    Polynomial v(1);

    u.coeff[2] = Number::ONE();
    u.coeff[1] = u1d;
    u.coeff[0] = u0d;

    v.coeff[1] = v1d;
    v.coeff[0] = v0d;

    Mumford ret(f, h, u, v);
    return ret;
}

Mumford Mumford::LangeAdd(const Mumford& m) const{
    std::cout << "Harley Addition." << std::endl;
    Polynomial u1 = this->u;
    Polynomial v1 = this->v;
    Polynomial u2 = m.u;
    Polynomial v2 = m.v;

    Number u11 = u1.coeff[1];
    Number u10 = u1.coeff[0];
    Number u21 = u2.coeff[1];
    Number u20 = u2.coeff[0];

    Number v11 = v1.coeff[1];
    Number v10 = v1.coeff[0];
    Number v21 = v2.coeff[1];
    Number v20 = v2.coeff[0];

    Number f6 = f.coeff[6];
    Number f5 = f.coeff[5];
    Number f4 = f.coeff[4];

    Polynomial u(2);
    Polynomial v(1);

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
        std::cout << "Special case." << std::endl;
        Mumford ret(f, h, u, v);
        return ret;
    }

    // 4. l' を計算．
    Number l3d = s1d;
    Number l2d = u21 * s1d + s0d;
    Number l1d = u21 * s0d + u20 * s1d;
    Number l0d = u20 * s0d;

    // 5. u' を計算．
    Number k4 = f6;
    Number k3 = f5 - f6 * u21;
    Number k2 = f4 - f6 * u20 - (f5 - f6 * u21) * u21;

    Number t4 = s1d * l3d - k4 * r * r;
    Number t3 = s1d * l2d + s0d * l3d - k3 * r * r;
    Number t2 = s1d * (l1d + r * v21 * 2) + s0d * l2d - k2 * r * r;

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
    Number l3 = l3d * w3;
    Number l2 = l2d * w3;
    Number l1 = l1d * w3;
    Number l0 = l0d * w3;
    Number v1d = -l1 - v21 + (u0d - u1d * u1d) * l3 + u1d * l2;
    Number v0d = -l0 - v20 - u1d * u0d * l3 + u0d * l2;

    u.coeff[2] = Number::ONE();
    u.coeff[1] = u1d;
    u.coeff[0] = u0d;

    v.coeff[1] = v1d;
    v.coeff[0] = v0d;

    Mumford ret(f, h, u, v);

    return ret;
}

Mumford Mumford::doubling(){
    Polynomial u = this->u;
    Polynomial v = this->v;

    Number u1 = u.coeff[1];
    Number u0 = u.coeff[0];
    Number v1 = v.coeff[1];
    Number v0 = v.coeff[0];

    Number h0 = this->h.coeff[0];
    Number h1 = this->h.coeff[1];
    Number h2 = this->h.coeff[2];

    Number f2 = this->f.coeff[2];
    Number f3 = this->f.coeff[3];
    Number f4 = this->f.coeff[4];

    // 1. v~ を計算．
    Number v1t = h1 + v1 * 2 - h2 * u1;
    Number v0t = h0 + v0 * 2 - h2 * u0;

    // 2. v~ と u の終結式を計算．
    Number w0 = v1 * v1;
    Number w1 = u1 * u1;
    Number w2 = v1t * v1t;
    Number w3 = u1 * v1t;
    Number r = u0 * w2 + v0t * (v0t - w3);

    // 3. r の almost inverse を計算．
    Number inv1d = v1t * (-1);
    Number inv0d = v0t - w3;

    // 4. k' を計算．
    w3 = f3 + w1;
    Number w4 = u0 * 2;
    Number k1d = (w1 - f4 * u1) * 2 + w3 - w4 - h2 * v1;
    Number k0d = u1 * (w4 * 2 - w3 + f4 * u1 + h2 * v1) + f2 - w0 - f4 * u0 * 2 - h1 * v1 - h2 * v0;

    // 5. s' を計算．
    w0 = k0d * inv0d;
    w1 = k1d * inv1d;
    Number s1d = (inv0d + inv1d) * (k0d + k1d) - w0 - w1 * (Number::ONE() + u1);
    Number s0d = w0 - u0 * w1;

    if(!s1d.isZero()){
        // 6. s'' を計算．
        w1 = (r * s1d).inv();
        w2 = r * w1;
        w3 = s1d * s1d * w1;
        w4 = r * w2;
        Number w5 = w4 * w4;
        Number s0dd = s0d * w2;

        // 7. l' を計算．
        Number l2d = u1 + s0dd;
        Number l1d = u1 * s0dd + u0;
        Number l0d = u0 * s0dd;

        // 8. u' を計算．
        Number u0d = s0dd * s0dd + w4 * (h2 * (s0dd - u1) + v1 * 2 + h1) + w5 * (u1 * 2 - f4);
        Number u1d = s0dd * 2 + h2 * w4 - w5;

        // 9. v' を計算．
        w1 = l2d - u1d;
        w2 = u1d * w1 + u0d - l1d;
        Number v1d = w2 * w3 - v1 - h1 + h2 * u1d;
        w2 = u0d * w1 - l0d;
        Number v0d = w2 * w3 - v0 - h0 + h2 * u0d;

        Polynomial u(2);
        Polynomial v(1);

        u.coeff[2] = Number::ONE();
        u.coeff[1] = u1d;
        u.coeff[0] = u0d;

        v.coeff[1] = v1d;
        v.coeff[0] = v0d;

        Mumford ret(f, h, u, v);
        return ret;
    }else{
        std::cout << "Special case." << std::endl;
        // サブルーチン

        // 6'. s を計算．
        w1 = r.inv();
        Number s0 = s0d * w1;
        w2 = u0 * s0 + v0 + h0;

        // 7'. u' を計算．
        Number u0d = f4 - s0 * s0 - s0 * h2 - u1 * 2;

        // 8. v' を計算．
        Number w1 = s0 * (u1 - u0d) - h2 * h2 * u0d + v1 + h1;
        Number v0d = v0d * w1 - w2;


        Polynomial u(1);
        Polynomial v(0);

        u.coeff[1] = Number::ONE();
        u.coeff[0] = u0d;

        v.coeff[0] = v0d;

        Mumford ret(f, h, u, v);
        return ret;
    }

    Mumford ret;
    return ret;
}

Mumford Mumford::CantorAdd(const Mumford& m) const{
    Polynomial h = this->h;

    Polynomial u1 = this->u;
    Polynomial v1 = this->v;
    Polynomial u2 = m.u;
    Polynomial v2 = m.v;

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

    Mumford ret(f, h, u, v);

    return ret;
}

Mumford Mumford::inv(){
    Polynomial f = this->f;
    Polynomial h = this->h;
    Polynomial u = this->u;
    Polynomial v = this->v;

    Mumford inv(f, h, u, (-(h + v)) % u);
    return inv;
}

Mumford Mumford::zero(){
    Polynomial f = this->f;
    Polynomial h = this->h;
    Polynomial ONE(0, 1);
    Polynomial ZERO(0, 0);

    Mumford zero(f, h, ONE, ZERO);
    return zero;
}

void Mumford::print(){
    std::cout << "[" << this->u << ", " << this->v << "]" << std::endl;
    return;
}
