#include "iostream"
#include "number.hpp"
#include "polynomial.hpp"
#include "mumford.hpp"

int main(){
    Number f0(-1);
    Number f1(3);
    Number f2(6);
    Number f3(-2);
    Number f4(-3);
    Number f5(1);

    Number h0(0);
    Number h1(1);
    Number h2(1);

    Number u10(5);
    Number u11(25);
    Number u12(1);

    Number v10(1);
    Number v11(4);

    Number u20(18);
    Number u21(13);
    Number u22(1);

    Number v20(-3);
    Number v21(-14);

    Polynomial f(5);
    Polynomial h(2);
    Polynomial u1(2);
    Polynomial v1(1);
    Polynomial u2(2);
    Polynomial v2(1);

    f.coeff[0] = f0; f.coeff[1] = f1; f.coeff[2] = f2; f.coeff[3] = f3; f.coeff[4] = f4; f.coeff[5] = f5;
    h.coeff[0] = h0; h.coeff[1] = h1; h.coeff[2] = h2;
    u1.coeff[0] = u10; u1.coeff[1] = u11; u1.coeff[2] = u12;
    v1.coeff[0] = v10; v1.coeff[1] = v11;
    u2.coeff[0] = u20; u2.coeff[1] = u21; u2.coeff[2] = u22;
    v2.coeff[0] = v20; v2.coeff[1] = v21;

    Mumford D1(f, h, u1, v1);
    Mumford D2(f, h, u2, v2);

    Mumford D = D1 + D2;
    D.print();

    /*
    auto tup = Polynomial::extended_gcd(u1, u2);
    Polynomial r = std::get<0>(tup);
    r.print();
    */
    
    return 0;
}
