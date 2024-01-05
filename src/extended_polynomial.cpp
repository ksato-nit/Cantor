#include "extended_polynomial.hpp"

ExtendedPolynomial::ExtendedPolynomial(){
    this->deg = 0;
    this->coeff.resize(1);
    this->coeff[0] = ExtendedNumber::ZERO();
}

ExtendedPolynomial::ExtendedPolynomial(int deg){
    this->deg = deg;
    this->coeff.resize(deg + 1);
}

ExtendedPolynomial::ExtendedPolynomial(int deg, int* coeff){
    this->deg = deg;
    this->coeff.resize(deg + 1);
    for(int i = 0; i <= deg; ++i){
        ExtendedNumber c(coeff[i], 0);
        this->coeff[i] = c;
    }
}

ExtendedPolynomial::ExtendedPolynomial(std::vector<ExtendedNumber> coeff){
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
    ExtendedNumber c(c0, 0);
    this->coeff[0] = c;
}

ExtendedPolynomial::ExtendedPolynomial(int deg, ExtendedNumber c0){
    this->deg = deg;
    this->coeff.resize(deg + 1);
    this->coeff[0] = c0;
}

ExtendedPolynomial::ExtendedPolynomial(int deg, ExtendedNumber c0, ExtendedNumber c1){
    this->deg = deg;
    this->coeff.resize(deg + 1);
    this->coeff[0] = c0;
    this->coeff[1] = c1;
}

bool ExtendedPolynomial::isZero() const{
    for(auto n : this->coeff){
        if(!n.isZero()){
            return false;
        }
    }
    return true;
}
