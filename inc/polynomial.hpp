#pragma once
#include "number.hpp"
#include "vector"
#include "iostream"

class Polynomial {
public:
    int deg;
    std::vector<Number> coeff;

    Polynomial(int deg);
    Polynomial(int deg, Number c0);
    Polynomial(int deg, Number c0, Number c1);
    Polynomial operator + (Polynomial q);
    Polynomial operator - (Polynomial q);
    Polynomial operator * (Polynomial q);
    void print();
};
