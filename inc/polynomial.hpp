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

    Number eval(Number x);

    Polynomial operator + (const Polynomial&) const;
    Polynomial operator - (const Polynomial&) const;
    Polynomial operator * (const Polynomial&) const;
    Polynomial operator * (const Number&) const;
    Polynomial operator / (const Polynomial& g) const;
    Polynomial operator % (const Polynomial&) const;
    Polynomial operator + () const;
    Polynomial operator - () const;
    bool operator == (const Polynomial& g) const;
    bool operator != (const Polynomial& g) const;
    static std::tuple<Polynomial, Polynomial> divide(Polynomial f, Polynomial g);
    static std::tuple<Polynomial, Polynomial, Polynomial> extended_gcd(Polynomial f, Polynomial g);
    void print();
};

inline std::ostream& operator << (std::ostream& os, const Polynomial& f){
    for(int i = f.deg; i >= 0; --i){
        os << f.coeff[i];
        if(i != 0){
            os << "x^" << i << " + ";
        }
    }
    return os;
}
