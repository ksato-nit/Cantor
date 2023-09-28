#include "mumford.hpp"

Mumford::Mumford(){
    this->u = Polynomial();
    this->v = Polynomial();
}

Mumford::Mumford(Polynomial u, Polynomial v){
    this->u = u;
    this->v = v;
}

Mumford Mumford::operator + (Mumford m){
    Mumford d;

    Polynomial u1 = this->u;
    Polynomial v1 = this->v;
    Polynomial u2 = m.u;
    Polynomial v2 = m.v;

    auto tup1 = Polynomial::extended_gcd(u1, u2);

    return d;
}

void Mumford::print(){
    this->u.print();
    this->v.print();
    return;
}
