#include "iostream"
#include "number.hpp"
#include "polynomial.hpp"

int main(){
    Number c0(2);
    Number c1(-13);
    Number c2(15);
    Number c3(1);

    Number b0(2);
    Number b1(5);
    Number b2(-12);

    Polynomial f(2);
    Polynomial g(2);
    f.coeff[0] = c0; f.coeff[1] = c1; f.coeff[2] = c2; // f.coeff[3] = c3;
    g.coeff[0] = b0; g.coeff[1] = b1; g.coeff[2] = b2;

    auto tup = Polynomial::extended_gcd(f, g);
    Polynomial r = std::get<0>(tup);
    r.print();

    return 0;
}
