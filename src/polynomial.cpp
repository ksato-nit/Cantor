#include "polynomial.hpp"

Polynomial::Polynomial(int deg){
    this->deg = deg;
}

Polynomial::Polynomial(int deg, Number c0){
    this->deg = deg;
    this->coeff.push_back(c0);
}

Polynomial::Polynomial(int deg, Number c0, Number c1){
    this->deg = deg;
    this->coeff.push_back(c0);
    this->coeff.push_back(c1);
}

Polynomial Polynomial::operator + (Polynomial q){
    // TODO : もう少しスマートな実装を考える．
    int deg = std::max(this->deg, q.deg);
    Polynomial r(deg);
    for(int i = 0; i <= deg; ++i){
        Number a = (this->deg >= i) ? this->coeff[i] : 0;
        Number b = (q.deg >= i) ? q.coeff[i] : 0;
        r.coeff.push_back(a + b);
    }
    return r;
}

void Polynomial::print(){
    for(int i = 0; i <= this->deg; ++i){
        std::cout << this->coeff[i].value << "x^" << i << ((i == this->deg) ? "" : " + ");
    }
    std::cout << std::endl;
}
