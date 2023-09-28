#pragma once
#include "number.hpp"
#include "polynomial.hpp"

class Mumford {
public:
    Polynomial u;
    Polynomial v;

    Polynomial f;
    Polynomial h;

    Mumford();
    Mumford(Polynomial f, Polynomial h);
    Mumford(Polynomial f, Polynomial h, Polynomial u, Polynomial v);

    Mumford operator + (Mumford q);
    void print();
};
