#include "polynomial.hpp"

Polynomial::Polynomial(int deg){
    this->deg = deg;
    this->coeff.resize(deg + 1);
}

Polynomial::Polynomial(int deg, Number c0){
    this->deg = deg;
    this->coeff.resize(deg + 1);
    this->coeff[0] = c0;
}

Polynomial::Polynomial(int deg, Number c0, Number c1){
    this->deg = deg;
    this->coeff.resize(deg + 1);
    this->coeff[0] = c0;
    this->coeff[1] = c1;
}

Polynomial Polynomial::operator + (Polynomial q){
    // TODO : もう少しスマートな実装を考える．
    int deg = std::max(this->deg, q.deg);
    Polynomial r(deg);
    for(int i = 0; i <= deg; ++i){
        Number a = (this->deg >= i) ? this->coeff[i] : 0;
        Number b = (q.deg >= i) ? q.coeff[i] : 0;
        r.coeff[i] = a + b;
    }
    return r;
}

Polynomial Polynomial::operator - (Polynomial q){
    int deg = std::max(this->deg, q.deg);
    Polynomial r(deg);
    for(int i = 0; i <= deg; ++i){
        Number a = (this->deg >= i) ? this->coeff[i] : 0;
        Number b = (q.deg >= i) ? q.coeff[i] : 0;
        r.coeff[i] = a - b;
    }
    return r;
}

Polynomial Polynomial::operator * (Polynomial q){
    int deg = this->deg + q.deg;
    Polynomial r(deg);
    for(int i = 0; i <= this->deg; ++i){
        for(int j = 0; j <= q.deg; ++j){
            Number a = (this->deg >= i) ? this->coeff[i] : 0;
            Number b = (q.deg >= j) ? q.coeff[j] : 0;
            r.coeff[i + j] = r.coeff[i + j] + a * b;
        }
    }
    return r;
}

void Polynomial::print(){
    for(int i = 0; i <= this->deg; ++i){
        std::cout << this->coeff[i].value << "x^" << i << ((i == this->deg) ? "" : " + ");
    }
    std::cout << std::endl;
}
