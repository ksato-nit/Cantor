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
    Mumford ret = this->CantorAdd(m);
    return ret;
}

Mumford Mumford::HarleyAdd(Mumford m){
    

    Mumford ret(f, h, u, v);
    return ret;
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
