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

void Polynomial::print(){
    for(int i = 0; i <= this->deg; ++i){
        std::cout << this->coeff[i].value << "x^" << i << ((i == this->deg) ? "" : " + ");
    }
    std::cout << std::endl;
}
