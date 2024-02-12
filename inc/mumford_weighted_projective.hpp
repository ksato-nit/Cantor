#pragma once
#include "number.hpp"
#include "extended_number.hpp"
#include "polynomial.hpp"
#include "divisor.hpp"

template <class T>
class WeightedProjectiveMumford {
    static const int GENUS = 2;
public:
    T U1;
    T U0;
    T V1;
    T V0;
    T Z1;
    T Z2;

    Polynomial<T> f;
    Polynomial<T> h;

    WeightedProjectiveMumford<T>();
    WeightedProjectiveMumford<T>(Polynomial<T> f, Polynomial<T> h);
    WeightedProjectiveMumford<T>(Polynomial<T> f, Polynomial<T> h, T, T, T, T, T, T);
    WeightedProjectiveMumford<T>(Polynomial<T> f, Polynomial<T> h, T, T, T, T);

    WeightedProjectiveMumford<T> operator + (const WeightedProjectiveMumford<T>& q) const;
    WeightedProjectiveMumford<T> CostelloAdd(const WeightedProjectiveMumford<T>& q) const;
    WeightedProjectiveMumford<T> inv();
    WeightedProjectiveMumford<T> zero();
    void print();
};
