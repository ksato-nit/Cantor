#include "extended_polynomial.hpp"

ExtendedPolynomial::ExtendedPolynomial(){
    this->deg = 0;
    this->coeff.resize(1);
    this->coeff[0] = Number::ZERO();
}

ExtendedPolynomial::ExtendedPolynomial(int deg){
    this->deg = deg;
    this->coeff.resize(deg + 1);
}

ExtendedPolynomial::ExtendedPolynomial(int deg, int* coeff){
    this->deg = deg;
    this->coeff.resize(deg + 1);
    for(int i = 0; i <= deg; ++i){
        Number c(coeff[i]);
        this->coeff[i] = c;
    }
}

ExtendedPolynomial::ExtendedPolynomial(std::vector<Number> coeff){
    int deg = coeff.size() + 1;
    this->deg = deg;
    this->coeff.resize(deg + 1);
    for(int i = 0; i <= deg; ++i){
        this->coeff[i] = coeff[i];
    }
}

ExtendedPolynomial::ExtendedPolynomial(int deg, int c0){
    this->deg = deg;
    this->coeff.resize(deg + 1);
    Number c(c0);
    this->coeff[0] = c;
}

ExtendedPolynomial::ExtendedPolynomial(int deg, Number c0){
    this->deg = deg;
    this->coeff.resize(deg + 1);
    this->coeff[0] = c0;
}

ExtendedPolynomial::ExtendedPolynomial(int deg, Number c0, Number c1){
    this->deg = deg;
    this->coeff.resize(deg + 1);
    this->coeff[0] = c0;
    this->coeff[1] = c1;
}

