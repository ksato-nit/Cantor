#pragma once
#include "extended_number.hpp"
#include "vector"
#include "iostream"
#include "tuple"

class ExtendedPolynomial {
public:
    int deg;
    std::vector<ExtendedNumber> coeff;

    ExtendedPolynomial();
    ExtendedPolynomial(int deg);
    ExtendedPolynomial(int deg, int c0);
    // TODO : テンプレートにしたい……．
    ExtendedPolynomial(int deg, int* coeff);
    ExtendedPolynomial(std::vector<ExtendedNumber> coeff);
    ExtendedPolynomial(int deg, ExtendedNumber c0);
    ExtendedPolynomial(int deg, ExtendedNumber c0, ExtendedNumber c1);


};
