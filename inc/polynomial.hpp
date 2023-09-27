#pragma once
#include "number.hpp"
#include "vector"
#include "iostream"

class Polynomial {
private:
    int deg;
    std::vector<Number> coeff;
public:
    Polynomial(int deg);
    Polynomial(int deg, Number c0);
    Polynomial(int deg, Number c0, Number c1);
    Polynomial operator + (Polynomial q);
    Polynomial operator - (Polynomial q);
    Polynomial operator * (Polynomial q);
    void print();
};
