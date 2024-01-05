#pragma once
#include "extended_number.hpp"
#include "extended_polynomial.hpp"
#include "divisor.hpp"

class ExtendedMumford {
    static const int GENUS = 2;
public:
    ExtendedNumber u2, u1, u0;
    ExtendedNumber v1, v0;
    ExtendedNumber U1, U0;

    ExtendedPolynomial f;
    ExtendedPolynomial h;
    ExtendedMumford(ExtendedPolynomial f, ExtendedPolynomial h);
    ExtendedMumford(ExtendedPolynomial f, ExtendedPolynomial h, ExtendedNumber u1, ExtendedNumber u0, ExtendedNumber v1, ExtendedNumber v0);
    ExtendedMumford(ExtendedPolynomial f, ExtendedPolynomial h, ExtendedNumber u1, ExtendedNumber u0, ExtendedNumber v1, ExtendedNumber v0, ExtendedNumber U1, ExtendedNumber U0);
    
    ExtendedMumford CostelloAdd(const ExtendedMumford& q) const;
    ExtendedMumford LangeAdd(const ExtendedMumford& q) const;
    
    ExtendedMumford CostelloDoubling() const;
    ExtendedMumford LangeDoubling() const;
};
