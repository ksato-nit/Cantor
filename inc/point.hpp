#pragma once
#include "vector"
#include "number.hpp"

class Point {
public:
    Number x, y;

    Point();
    Point(Number, Number);
    Point involution();
};
