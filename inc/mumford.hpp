#pragma once
#include "number.hpp"
#include "polynomial.hpp"
#include "divisor.hpp"
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>

class Mumford {
    static const int GENUS = 2;
public:
    Number u2, u1, u0;
    Number v1, v0;
    Number U1, U0;
    //Polynomial u;
    //Polynomial v;

    Polynomial f;
    Polynomial h;

    Mumford();
    Mumford(Polynomial f, Polynomial h);
    Mumford(Polynomial f, Polynomial h, Number u1, Number u0, Number v1, Number v0);
    Mumford(Polynomial f, Polynomial h, Number u1, Number u0, Number v1, Number v0, Number U1, Number U0);
    Mumford(Polynomial f, Polynomial h, Divisor d);

    Mumford operator + (const Mumford& q) const;

    Mumford operator * (const boost::multiprecision::int1024_t& k) const;

    Mumford CostelloScalarMultiple(const boost::multiprecision::int1024_t& k) const;

    Mumford CantorAdd(const Mumford& q) const;
    
    Mumford CostelloAdd(const Mumford& q) const;
    Mumford LangeAdd(const Mumford& q) const;

    // 5 次・モニックにのみ対応
    Mumford HarleyAdd(const Mumford& q) const;
    Mumford HarleyAddDegenerated(const Mumford& q) const;
    
    Mumford CostelloDoubling() const;
    Mumford LangeDoubling() const;

    Mumford inv();
    Mumford zero() const; // todo: static で書き直したい
    static Mumford zero(const Polynomial& f, const Polynomial& h);
    void print() const;
    bool isZero() const;
};
