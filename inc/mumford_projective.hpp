#pragma once
#include "number.hpp"
#include "polynomial.hpp"
#include "divisor.hpp"

template <class T>
class ProjectiveMumford {
    static const int GENUS = 2;
public:
    T U1;
    T U0;
    T V1;
    T V0;
    T Z;
    T W1;
    T W0;

    Polynomial<T> f;
    Polynomial<T> h;

    ProjectiveMumford<T>();
    ProjectiveMumford<T>(Polynomial<T> f, Polynomial<T> h);
    ProjectiveMumford<T>(Polynomial<T> f, Polynomial<T> h, T, T, T, T, T, T, T);
    ProjectiveMumford<T>(Polynomial<T> f, Polynomial<T> h, T, T, T, T, T);
    ProjectiveMumford<T>(Polynomial<T> f, Polynomial<T> h, T, T, T, T);

    ProjectiveMumford<T> operator + (const ProjectiveMumford<T>& q) const;
    ProjectiveMumford<T> operator * (const mpz_class& k) const;

    ProjectiveMumford<T> CostelloAdd(const ProjectiveMumford<T>& q) const;
    ProjectiveMumford<T> LangeAdd(const ProjectiveMumford<T>& q) const;
    ProjectiveMumford<T> LangeDoubling() const;
    ProjectiveMumford<T> inv();
    ProjectiveMumford<T> zero() const;
    static ProjectiveMumford<T> zero(const Polynomial<T>& f, const Polynomial<T>& h);
    void print() const;
    bool isZero() const;
};
