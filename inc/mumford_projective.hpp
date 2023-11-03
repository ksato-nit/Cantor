#pragma once
#include "number.hpp"
#include "polynomial.hpp"
#include "divisor.hpp"

class ProjectiveMumford {
    static const int GENUS = 2;
public:
    Number U1;
    Number U0;
    Number V1;
    Number V0;
    Number Z;

    Polynomial f;
    Polynomial h;

    ProjectiveMumford();
    ProjectiveMumford(Polynomial f, Polynomial h);
    ProjectiveMumford(Polynomial f, Polynomial h, Number, Number, Number, Number, Number);
    ProjectiveMumford(Polynomial f, Polynomial h, Number, Number, Number, Number);

    ProjectiveMumford operator + (const ProjectiveMumford& q) const;
    ProjectiveMumford CostelloAdd(const ProjectiveMumford& q) const;
    ProjectiveMumford inv();
    ProjectiveMumford zero();
    void print();
};
