#pragma once
#include "number.hpp"
#include "polynomial.hpp"

class Mumford {
public:
    Polynomial u;
    Polynomial v;

    Mumford();
    Mumford(Polynomial u, Polynomial v);

    Mumford operator + (Mumford q);
    void print();
};
