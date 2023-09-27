#include "iostream"
#include "number.hpp"
#include "polynomial.hpp"

int main(){
    int deg = 1;
    Number c0(7);
    Number c1(3);

    Polynomial f(deg, c0, c1);

    f.print();

    return 0;
}