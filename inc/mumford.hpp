#pragma once
#include "number.hpp"
#include "polynomial.hpp"
#include "divisor.hpp"

class Mumford {
    static const int GENUS = 2;
public:
    Polynomial u;
    Polynomial v;

    Polynomial f;
    Polynomial h;

    Mumford();
    Mumford(Polynomial f, Polynomial h);
    Mumford(Polynomial f, Polynomial h, Polynomial u, Polynomial v);
    Mumford(Polynomial f, Polynomial h, Divisor d);

    Mumford operator + (const Mumford& q) const;

    Mumford operator * (const Number& k) const;

    Mumford CantorAdd(const Mumford& q) const;
    
    Mumford CostelloAdd(const Mumford& q) const;
    Mumford LangeAdd(const Mumford& q) const;

    // 5 次・モニックにのみ対応
    Mumford HarleyAdd(const Mumford& q) const;
    Mumford HarleyAddDegenerated(const Mumford& q) const;
    
    Mumford LangeDoubling() const;

    Mumford inv();
    Mumford zero() const; // todo: static で書き直したい
    void print();
};
