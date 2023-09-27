#include "iostream"
#include "number.hpp"
#include "polynomial.hpp"

int main(){
    int deg = 1;
    Number c0(7);
    Number c1(3);
    Number b0(5);
    Number b1(2);

    Polynomial f(deg, c0, c1);
    Polynomial g(deg, b0, b1);

    Polynomial h = f + g;

    h.print();

    return 0;
}