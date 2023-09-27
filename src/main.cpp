#include "iostream"
#include "number.hpp"
#include "polynomial.hpp"

int main(){
    Number c0(7);
    Number c1(-2);
    Number c2(4);
    Number c3(1);

    Number b0(2);
    Number b1(-3);
    Number b2(1);

    Polynomial f(3);
    Polynomial g(2);
    f.coeff[0] = c0; f.coeff[1] = c1; f.coeff[2] = c2; f.coeff[3] = c3;
    g.coeff[0] = b0; g.coeff[1] = b1; g.coeff[2] = b2;

    std::tuple<Polynomial, Polynomial> tup = Polynomial::divide(f, g);

    return 0;
}
