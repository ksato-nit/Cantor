#pragma once
#include "number.hpp"
#include "vector"
#include "iostream"
#include "tuple"

class Polynomial {
public:
    int deg;
    std::vector<Number> coeff;

    Polynomial();
    Polynomial(int deg);
    Polynomial(int deg, int c0);
    // TODO : テンプレートにしたい……．
    Polynomial(int deg, int* coeff);
    Polynomial(std::vector<Number> coeff);
    Polynomial(int deg, Number c0);
    Polynomial(int deg, Number c0, Number c1);

    void resize(int deg);
    void resize(std::vector<Number> coeff);

    bool isZero();

    void normalize();

    Polynomial operator + (Polynomial q);
    Polynomial operator - (Polynomial q);
    Polynomial operator * (Polynomial q);
    Polynomial operator * (Number k);
    Polynomial operator / (Polynomial g);
    Polynomial operator % (Polynomial g);
    static std::tuple<Polynomial, Polynomial> divide(Polynomial f, Polynomial g);
    static std::tuple<Polynomial, Polynomial, Polynomial> extended_gcd(Polynomial f, Polynomial g);
    void print();
};
