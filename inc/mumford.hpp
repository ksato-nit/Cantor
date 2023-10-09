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

    Mumford operator + (Mumford q);
    Mumford CantorAdd(Mumford q);
    Mumford HarleyAdd(Mumford q);
    Mumford inv();
    Mumford zero();
    void print();
};
