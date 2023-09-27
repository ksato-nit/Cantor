#pragma once
#include "number.hpp"
#include "vector"
#include "iostream"
#include "tuple"

class Polynomial {
public:
    int deg;
    std::vector<Number> coeff;

    Polynomial(int deg);
    Polynomial(int deg, Number c0);
    Polynomial(int deg, Number c0, Number c1);
    // TODO : std::vector を渡して初期化できるようにする．

    Polynomial operator + (Polynomial q);
    Polynomial operator - (Polynomial q);
    Polynomial operator * (Polynomial q);
    Polynomial operator * (Number k);
    static std::tuple<Polynomial, Polynomial> divide(Polynomial f, Polynomial g);
    static std::tuple<Polynomial, Polynomial, Polynomial> extended_gcd(Polynomial f, Polynomial g);
    void print();
};
