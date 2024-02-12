#pragma once
#include "number.hpp"
#include "vector"
#include "iostream"
#include "tuple"

template <class T>
class Polynomial {
public:
    int deg;
    std::vector<T> coeff;

    Polynomial<T>();
    Polynomial<T>(int deg);
    Polynomial<T>(int deg, int c0);
    Polynomial<T>(int deg, int* coeff);
    Polynomial<T>(std::vector<T> coeff);
    Polynomial<T>(int deg, T c0);
    Polynomial<T>(int deg, T c0, T c1);

    void resize(int deg);
    void resize(std::vector<T> coeff);

    bool isZero();

    void normalize();

    T eval(T x) const;

    Polynomial<T> operator + (const Polynomial<T>&) const;
    Polynomial<T> operator - (const Polynomial<T>&) const;
    Polynomial<T> operator * (const Polynomial<T>&) const;
    Polynomial<T> operator * (const T&) const;
    Polynomial<T> operator / (const Polynomial<T>& g) const;
    Polynomial<T> operator % (const Polynomial<T>&) const;
    Polynomial<T> operator + () const;
    Polynomial<T> operator - () const;
    bool operator == (const Polynomial<T>& g) const;
    bool operator != (const Polynomial<T>& g) const;
    static std::tuple<Polynomial<T>, Polynomial<T>> divide(Polynomial<T> f, Polynomial<T> g);
    static std::tuple<Polynomial<T>, Polynomial<T>, Polynomial<T>> extended_gcd(Polynomial<T> f, Polynomial<T> g);
    void print() const;

    Polynomial<T> derivative() const;
};

template <typename T>
inline std::ostream& operator << (std::ostream& os, const Polynomial<T>& f){
    for(int i = f.deg; i >= 0; --i){
        os << f.coeff[i];
        if(i != 0){
            os << "x^" << i << " + ";
        }
    }
    return os;
}
