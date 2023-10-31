#include "mumford_projective.hpp"

ProjectiveMumford::ProjectiveMumford(){
    this->f = Polynomial();
    this->h = Polynomial();
    this->U1 = Number::ZERO();
    this->U0 = Number::ONE();
    this->V1 = Number::ZERO();
    this->V0 = Number::ZERO();
    this->Z = Number::ZERO();
}

ProjectiveMumford::ProjectiveMumford(Polynomial f, Polynomial h){
    this->f = f;
    this->h = h;
    this->U1 = Number::ZERO();
    this->U0 = Number::ONE();
    this->V1 = Number::ZERO();
    this->V0 = Number::ZERO();
    this->Z = Number::ZERO();
}

ProjectiveMumford::ProjectiveMumford(Polynomial f, Polynomial h, Number U1, Number U0, Number V1, Number V0, Number Z){
    this->f = f;
    this->h = h;
    this->U1 = U1;
    this->U0 = U0;
    this->V1 = V1;
    this->V0 = V0;
    this->Z = Z;
}

ProjectiveMumford::ProjectiveMumford(Polynomial f, Polynomial h, Number U1, Number U0, Number V1, Number V0){
    this->f = f;
    this->h = h;
    this->U1 = U1;
    this->U0 = U0;
    this->V1 = V1;
    this->V0 = V0;
    this->Z = Number::ONE();
}

ProjectiveMumford ProjectiveMumford::operator + (ProjectiveMumford m){
    // deg u1 = deg u2 = 2
    ProjectiveMumford ret = this->HarleyAdd(m);
    return ret;
}

ProjectiveMumford ProjectiveMumford::HarleyAdd(ProjectiveMumford m){
    std::cout << "Harley Addition." << std::endl;

    /*
    Number u11 = u1.coeff[1];
    Number u10 = u1.coeff[0];
    Number u21 = u2.coeff[1];
    Number u20 = u2.coeff[0];

    Number v11 = v1.coeff[1];
    Number v10 = v1.coeff[0];
    Number v21 = v2.coeff[1];
    Number v20 = v2.coeff[0];

    Number h0 = this->h.coeff[0];
    Number h1 = this->h.coeff[1];
    Number h2 = this->h.coeff[2];

    Number f4 = this->f.coeff[4];

    // 1. u1, u2 の終結式を計算．
    Number z1 = u11 - u21;
    Number z2 = u20 - u10;
    Number z3 = u11 * z1 + z2;
    Number r = z2 * z3 + z1 * z1 * u10;

    // 2. u2 の almost inverse (mod u1) を計算．
    Number inv1 = z1;
    Number inv0 = z3;

    // 3. s' を計算．
    Number w0 = v10 - v20;
    Number w1 = v11 - v21;
    Number w2 = inv0 * w0;
    Number w3 = inv1 * w1;
    Number s1d = (inv0 + inv1) * (w0 + w1) - w2 - w3 * (Number::ONE() + u11);
    Number s0d = w2 - u10 * w3;

    if(!s1d.isZero()){
        // 4. s'' を計算．
        Number w1 = (r * s1d).inv();
        Number w2 = r * w1;
        Number w3 = s1d * s1d * w1;
        Number w4 = r * w2;
        Number w5 = w4 * w4;
        Number s0dd = s0d * w2;

        // 5. l' を計算．
        Number l2d = u21 + s0dd;
        Number l1d = u21 * s0dd + u20;
        Number l0d = u20 * s0dd;


        // 6. u' を計算．
        Number u0d = (s0dd - u11) * (s0dd - z1 + h2 * w4) - u10 + l1d + (h1 + v21 * 2) * w4 + (u21 * 2 + z1 - f4) * w5;
        Number u1d = s0dd * 2 - z1 + h2 * w4 - w5;

        // 7. v' を計算．
        w1 = l2d - u1d;
        w2 = u1d * w1 + u0d - l1d;
        Number v1d = w2 * w3 - v21 - h1 + h2 * u1d;
        w2 = u0d * w1 - l0d;
        Number v0d = w2 * w3 - v20 - h0 + h2 * u0d;

        Polynomial u(2);
        Polynomial v(1);

        u.coeff[2] = Number::ONE();
        u.coeff[1] = u1d;
        u.coeff[0] = u0d;

        v.coeff[1] = v1d;
        v.coeff[0] = v0d;

        Mumford ret(f, h, u, v);
        return ret;
    }else{
        std::cout << "Special case." << std::endl;
        // サブルーチン

        // 4'. s を計算．
        Number inv = r.inv();
        Number s0 = s0d * inv;

        // 5'. u' を計算．
        Number u0d = f4 - u21 - u11 - s0*s0 - s0 * h2;

        // 6'. v' を計算．
        Number w1 = s0 * (u21 + u0d) + h1 + v21 + h2 * u0d;
        Number w2 = s0 + v20 + h0;
        Number v0d = u0d * w1 - w2;

        Polynomial u(1);
        Polynomial v(0);

        u.coeff[1] = Number::ONE();
        u.coeff[0] = u0d;

        v.coeff[0] = v0d;

        Mumford ret(f, h, u, v);
        return ret;
    }
    */

    ProjectiveMumford ret(f, h);
    return ret;
}

ProjectiveMumford ProjectiveMumford::inv(){
    Polynomial f = this->f;
    Polynomial h = this->h;
    Number U1 = this->U1;
    Number U0 = this->U0;
    Number V1 = this->V1;
    Number V0 = this->V0;
    Number Z = this->Z;

    Number h2 = this->h.coeff[2];
    Number h1 = this->h.coeff[1];
    Number h0 = this->h.coeff[0];
    // - (h + v) % u を計算．
    Number V1d = U1 * h2 - (V1 + h1) * Z;
    Number V0d = U0 * h2 - (V0 + h0) * Z;
    ProjectiveMumford inv(f, h, U1, U0, V1d, V0d, Z);
    return inv;
}

ProjectiveMumford ProjectiveMumford::zero(){
    Polynomial f = this->f;
    Polynomial h = this->h;

    ProjectiveMumford zero(f, h);
    return zero;
}

void ProjectiveMumford::print(){
    std::cout << "[";
    std::cout << this->U1.value << " " << this->U0.value;
    std::cout << ", ";
    std::cout << this->V1.value << " " << this->V0.value;
    std::cout << ", ";
    std::cout << this->Z.value;
    std::cout << "]" << std::endl;
    return;
}
