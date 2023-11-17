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
    Mumford CantorAdd(const Mumford& q) const;
    Mumford HarleyAdd(const Mumford& q) const;
    Mumford CostelloAdd(const Mumford& q) const;
    Mumford HarleyAddDegenerated(const Mumford& q) const;
    Mumford LangeAdd(const Mumford& q) const;
    Mumford doubling();
    Mumford inv();
    Mumford zero();
    void print();
};
