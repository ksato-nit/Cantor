#pragma once
#include "number.hpp"

class Polynomial {
private:
    Number x, y;
public:
    Polynomial();
    Polynomial operator + (Polynomial q);
    Polynomial operator - (Polynomial q);
    Polynomial operator * (Polynomial q);
};
