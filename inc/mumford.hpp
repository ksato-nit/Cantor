#pragma once
#include "number.hpp"
#include "polynomial.hpp"
#include "divisor.hpp"

template <class T>
class Mumford {
    static const int GENUS = 2;
public:
    T u2, u1, u0;
    T v1, v0;

    Polynomial<T> f;
    Polynomial<T> h;

    Mumford();
    Mumford(Polynomial<T> f, Polynomial<T> h);
    Mumford(Polynomial<T> f, Polynomial<T> h, T u1, T u0, T v1, T v0);
    Mumford(Polynomial<T> f, Polynomial<T> h, Divisor d);

    Mumford<T> operator + (const Mumford<T>& q) const;

    Mumford<T> operator * (const mpz_class& k) const;

    Mumford<T> CostelloScalarMultiple(const mpz_class& k) const;

    Mumford<T> CantorAdd(const Mumford<T>& q) const;
    
    Mumford<T> CostelloAdd(const Mumford<T>& q) const;
    Mumford<T> LangeAdd(const Mumford<T>& q) const;

    // 5 次・モニックにのみ対応
    Mumford<T> HarleyAdd(const Mumford<T>& q) const;
    Mumford<T> HarleyAddDegenerated(const Mumford<T>& q) const;
    
    Mumford<T> CostelloDoubling() const;
    Mumford<T> LangeDoubling() const;

    Mumford<T> inv();
    Mumford<T> zero() const; // todo: static で書き直したい
    static Mumford<T> zero(const Polynomial<T>& f, const Polynomial<T>& h);
    void print() const;
    bool isZero() const;
};
