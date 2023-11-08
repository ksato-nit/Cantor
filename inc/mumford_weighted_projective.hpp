#pragma once
#include "number.hpp"
#include "polynomial.hpp"
#include "divisor.hpp"

class WeightedProjectiveMumford {
    static const int GENUS = 2;
public:
    Number U1;
    Number U0;
    Number V1;
    Number V0;
    Number Z1;
    Number Z2;

    Polynomial f;
    Polynomial h;

    WeightedProjectiveMumford();
    WeightedProjectiveMumford(Polynomial f, Polynomial h);
    WeightedProjectiveMumford(Polynomial f, Polynomial h, Number, Number, Number, Number, Number, Number);
    WeightedProjectiveMumford(Polynomial f, Polynomial h, Number, Number, Number, Number);

    WeightedProjectiveMumford operator + (const WeightedProjectiveMumford& q) const;
    WeightedProjectiveMumford CostelloAdd(const WeightedProjectiveMumford& q) const;
    WeightedProjectiveMumford inv();
    WeightedProjectiveMumford zero();
    void print();
};
