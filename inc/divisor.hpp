#pragma once
#include "vector"

#include "point.hpp"

class Divisor {
public:
    std::vector< std::pair<Point, int> > points;

    Divisor();
    Divisor operator + (Divisor d2);
    Divisor operator - (Divisor d2);
    Divisor inv();
};
