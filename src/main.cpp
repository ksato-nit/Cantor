#include "iostream"
#include "number.hpp"
#include "polynomial.hpp"
#include "mumford.hpp"

int main(){
    int fc[6] = {-1, 3, 6, -2, -3, 1};
    int hc[3] = {0, 1, 1};
    int u1c[3] = {5, 25, 1};
    int v1c[2] = {-30, -27};
    int u2c[3] = {28, 25, 1};
    int v2c[2] = {-26, -18};

    Polynomial f(5, fc);
    Polynomial h(2, hc);
    Polynomial u1(2, u1c);
    Polynomial v1(1, v1c);
    Polynomial u2(2, u2c);
    Polynomial v2(1, v2c);

    Mumford D1(f, h, u1, v1);
    Mumford D2(f, h, u2, v2);

    Mumford now = D1;
    for(int i = 1; i < 1012; ++i){
        Mumford D1_new = now + D1;
        D1_new.print();
        now = D1_new;
    }

    return 0;
}
