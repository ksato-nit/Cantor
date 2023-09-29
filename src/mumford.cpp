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

Mumford Mumford::operator + (Mumford m){
    Polynomial h = this->h;

    Polynomial u1 = this->u;
    Polynomial v1 = this->v;
    Polynomial u2 = m.u;
    Polynomial v2 = m.v;

    auto tup1 = Polynomial::extended_gcd(u1, u2);
    Polynomial d1 = std::get<0>(tup1);
    Polynomial e1 = std::get<1>(tup1);
    Polynomial e2 = std::get<2>(tup1);

    /*
    d1.print();
    e1.print();
    e2.print();

    puts("Step 2");

    Polynomial vv = v1 + v2 + h;
    vv.print();
    */

    auto tup2 = Polynomial::extended_gcd(d1, v1 + v2 + h);
    Polynomial d = std::get<0>(tup2);
    Polynomial c1 = std::get<1>(tup2);
    Polynomial c2 = std::get<2>(tup2);

    /*
    d.print();
    c1.print();
    c2.print();

    Polynomial dd = (c1 * d1) + (c2 * (v1 + v2 + h));
    dd.print();
    */

    Polynomial s1 = c1 * e1;
    Polynomial s2 = c1 * e2;

    Polynomial u = u1 * u2 / (d * d);
    Polynomial v = ((s1 * u1 * v2 + s2 * u2 * v1 + c2 * (v1 * v2 + f)) / d) % u;

    u.print();
    v.print();

    Number MINUS_ONE(-1);

    do{
        Polynomial ud = (f - v * h - (v * v)) / u;
        Polynomial vr = ((h + v) * MINUS_ONE) % ud;

        u = ud;
        v = vr;
    }while(u.deg < GENUS);

    u.print();
    v.print();

    //Number lc = u.coeff[u.deg];

    //u = u * lc.inv();

    Mumford ret(f, h, u, v);

    return ret;
}

void Mumford::print(){
    this->u.print();
    this->v.print();
    return;
}
