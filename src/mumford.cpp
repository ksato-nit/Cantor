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
            Polynomial u(1, Number::ONE(), p.first.x * Number::MINUS_ONE());
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

        std::vector<Number> uc = {(p1.x * p2.x), (p1.x + p2.x) * Number::MINUS_ONE(), Number::ONE()};
        Polynomial u(uc);
        this->u = u;

        Number c1 = (p1.y - p2.y) / (p1.x - p2.x);
        Number c0 = (p1.x * p2.y - p2.x * p1.y) / (p1.x - p2.x);
        Polynomial v(1, c0, c1);
        this->v = v;

    }
}

Mumford Mumford::operator + (Mumford m){
    Mumford ret = this->HarleyAdd(m);
    return ret;
}

Mumford Mumford::HarleyAdd(Mumford m){
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

Mumford Mumford::CantorAdd(Mumford m){
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
        Polynomial vr = ((h + v)*(Number::MINUS_ONE())) % ud;

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

    Mumford inv(f, h, u, ((h + v) * Number::MINUS_ONE()) % u);
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
    std::cout << "[";
    this->u.print();
    std::cout << ", ";
    this->v.print();
    std::cout << "]" << std::endl;
    return;
}
