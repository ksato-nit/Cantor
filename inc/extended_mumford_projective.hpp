#pragma once
#include "extended_number.hpp"
#include "extended_polynomial.hpp"

class ExtendedProjectiveMumford {
    static const int GENUS = 2;
public:
    ExtendedNumber U1;
    ExtendedNumber U0;
    ExtendedNumber V1;
    ExtendedNumber V0;
    ExtendedNumber Z;
    ExtendedNumber W1;
    ExtendedNumber W0;

    ExtendedPolynomial f;
    ExtendedPolynomial h;

    ExtendedProjectiveMumford();
    ExtendedProjectiveMumford(ExtendedPolynomial f, ExtendedPolynomial h);
    ExtendedProjectiveMumford(ExtendedPolynomial f, ExtendedPolynomial h, ExtendedNumber, ExtendedNumber, ExtendedNumber, ExtendedNumber, ExtendedNumber, ExtendedNumber, ExtendedNumber);
    ExtendedProjectiveMumford(ExtendedPolynomial f, ExtendedPolynomial h, ExtendedNumber, ExtendedNumber, ExtendedNumber, ExtendedNumber, ExtendedNumber);
    ExtendedProjectiveMumford(ExtendedPolynomial f, ExtendedPolynomial h, ExtendedNumber, ExtendedNumber, ExtendedNumber, ExtendedNumber);

    ExtendedProjectiveMumford operator * (const mpz_class& k) const;

    ExtendedProjectiveMumford LangeAdd(const ExtendedProjectiveMumford& q) const;
    ExtendedProjectiveMumford LangeDoubling() const;
    static ExtendedProjectiveMumford zero(const ExtendedPolynomial& f, const ExtendedPolynomial& h);
    void print() const;
    bool isZero() const;
};
