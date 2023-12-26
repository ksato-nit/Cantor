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
    Number W1;
    Number W0;

    Polynomial f;
    Polynomial h;

    ProjectiveMumford();
    ProjectiveMumford(Polynomial f, Polynomial h);
    ProjectiveMumford(Polynomial f, Polynomial h, Number, Number, Number, Number, Number, Number, Number);
    ProjectiveMumford(Polynomial f, Polynomial h, Number, Number, Number, Number, Number);
    ProjectiveMumford(Polynomial f, Polynomial h, Number, Number, Number, Number);

    ProjectiveMumford operator + (const ProjectiveMumford& q) const;
    ProjectiveMumford operator * (const NTL::ZZ& k) const;

    ProjectiveMumford CostelloAdd(const ProjectiveMumford& q) const;
    ProjectiveMumford LangeAdd(const ProjectiveMumford& q) const;
    ProjectiveMumford LangeDoubling() const;
    ProjectiveMumford inv();
    ProjectiveMumford zero() const;
    static ProjectiveMumford zero(const Polynomial& f, const Polynomial& h);
    void print() const;
    bool isZero() const;
};
