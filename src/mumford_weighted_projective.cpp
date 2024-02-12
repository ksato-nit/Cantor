#include "mumford_weighted_projective.hpp"

template <class T>
WeightedProjectiveMumford<T>::WeightedProjectiveMumford(){
    this->f = Polynomial<T>();
    this->h = Polynomial<T>();
    this->U1 = T::ZERO();
    this->U0 = T::ONE();
    this->V1 = T::ZERO();
    this->V0 = T::ZERO();
    this->Z1 = T::ONE();
    this->Z2 = T::ONE();
}

template <class T>
WeightedProjectiveMumford<T>::WeightedProjectiveMumford(Polynomial<T> f, Polynomial<T> h){
    this->f = f;
    this->h = h;
    this->U1 = T::ZERO();
    this->U0 = T::ONE();
    this->V1 = T::ZERO();
    this->Z1 = T::ONE();
    this->Z2 = T::ONE();
}

template <class T>
WeightedProjectiveMumford<T>::WeightedProjectiveMumford(Polynomial<T> f, Polynomial<T> h, T U1, T U0, T V1, T V0, T Z1, T Z2){
    this->f = f;
    this->h = h;
    this->U1 = U1;
    this->U0 = U0;
    this->V1 = V1;
    this->V0 = V0;
    this->Z1 = Z1;
    this->Z2 = Z2;
}

template <class T>
WeightedProjectiveMumford<T>::WeightedProjectiveMumford(Polynomial<T> f, Polynomial<T> h, T U1, T U0, T V1, T V0){
    this->f = f;
    this->h = h;
    this->U1 = U1;
    this->U0 = U0;
    this->V1 = V1;
    this->V0 = V0;
    this->Z1 = T::ONE();
    this->Z2 = T::ONE();
}

template <class T>
WeightedProjectiveMumford<T> WeightedProjectiveMumford<T>::operator + (const WeightedProjectiveMumford<T>& m) const{
    // deg u1 = deg u2 = 2
    WeightedProjectiveMumford<T> ret = this->CostelloAdd(m);
    return ret;
}

template <class T>
WeightedProjectiveMumford<T> WeightedProjectiveMumford<T>::CostelloAdd(const WeightedProjectiveMumford<T>& m) const{
    std::cout << "Weighted Projective Costello Addition." << std::endl;
    T U11 = this->U1;
    T U10 = this->U0;
    T U21 = m.U1;
    T U20 = m.U0;

    T V11 = this->V1;
    T V10 = this->V0;
    T V21 = m.V1;
    T V20 = m.V0;

    T Z11 = this->Z1;
    T Z12 = this->Z2;
    T Z21 = m.Z1;
    T Z22 = m.Z2;

    T f6 = this->f.coeff[6];
    T f5 = this->f.coeff[5];
    T f4 = this->f.coeff[4];

    // 17M, 8S
    T Z11S = Z11 * Z11;
    T Z21S = Z21 * Z21;
    T Z11Q = Z11S * Z11S;
    T Z21Q = Z21S * Z21S;    
    T Z12S = Z12 * Z12;
    T Z22S = Z22 * Z22;
    T Z11Z12 = Z11 * Z12;
    T ZALL = Z11Z12 * Z21 * Z22;
    T U11Z21S = U11 * Z21S;
    T U21Z11S = U21 * Z11S;
    T U10Z21S = U10 * Z21S;
    T U20Z11S = U20 * Z11S;
    T Z21TZ22 = Z21S * Z21 * Z22;
    T Z11TZ12 = Z11S * Z11 * Z12;
    T V11Z21TZ22 = V11 * Z21TZ22;
    T V21Z11TZ12 = V21 * Z11TZ12;
    T V10Z21TZ22 = V10 * Z21TZ22;
    T V20Z11TZ12 = V20 * Z11TZ12;

    T ZZAS = Z11S * Z21S;
    T ZZBS = Z12S * Z22S;

    T T11 = U11 * U11;
    T T21 = U21 * U21;

    // 10M

    // (Z11Z21)^4 がかかっている．
    T M1 = (T11 - U10 * Z11S) * Z21Q - (T21 - U20 * Z21S) * Z11Q;
    T M2 = U21Z11S * U20Z11S - U11Z21S * U10Z21S;
    // (Z11Z21)^3 Z12 Z22 がかかっている．
    T M3 = U11Z21S - U21Z11S;
    T M4 = U20Z11S - U10Z21S;

    T z1 = V10Z21TZ22 - V20Z11TZ12;
    T z2 = V11Z21TZ22 - V21Z11TZ12;

    // (Z11Z21)^7 Z12 Z22 がかかっている．
    T t1 = (M2 - z1) * (z2 - M1);
    T t2 = (-z1 - M2) * (z2 + M1);
    // (Z11Z21)^5 Z12 Z22 がかかっている．
    T t3 = (-z1 + M4) * (z2 - M3);
    T t4 = (-z1 - M4) * (z2 + M3);


    // 1M
    // (Z11Z21)^7 Z12 Z22 がかかっている．
    T l2_num = t1 - t2;
    // (Z11Z21)^5 Z12 Z22 がかかっている．
    T l3_num = t3 - t4;
    // (Z11Z21)^6 がかかっている．
    T d = -(M2 - M4) * (M1 + M3);
    d = d + d + (t3 + t4) - t1 - t2;

    // 3M, 2S
    T A = d * d;
    T B = ZZAS * l3_num * l3_num - f6 * A * ZZBS;

    // 47M
    T Ud1 = -B * (U11 * Z21S + U21 * Z11S) - (f5 * A * ZZBS - l2_num * l3_num - l2_num * l3_num) * ZZAS;
    T Zd1 = B * ZZAS;

    T Ud0 = l3_num * ZZAS * (l3_num * ZZAS * (U10 * Z11S - T11) * Z12 + (l2_num * U11 * Z11Z12 + V11 * ZALL * d) * Z11) * 2;
    Ud0 = Ud0 + (l2_num * l2_num - f4 * ZZAS * ZZBS * A) * Z11Q * Z12;
    Ud0 = Ud0 * Z21S * Zd1;
    Ud0 = Ud0 - ((U11 * U21 + (U10 * Z21S + U20 * Z11S)) * Zd1 + (U11 * Z21S + U21 * Z11S) * Ud1) * Z11S * Z12 * B * ZZAS;

    Ud0 = Ud0 * Z12;
    Ud1 = Ud1 * B * Z11Q * Z11S * Z21Q * Z12S;
    Zd1 = Z11Q * Z21S * Z21 * B * Z12;

    // 2S
    T Zd1S = Zd1 * Zd1;
    T Zd1Q = Zd1S * Zd1S;

    // 30M
    T c3 = (Ud0 * Zd1S - Ud1 * Ud1) * Z11Q - (U10 * Z11S - T11) * Zd1Q;
    T c2 = (Ud1 * Z11S - U11 * Zd1S) * (Z11S * Zd1S);
    T c1 = Z11 * Zd1Q * d * ZALL;
    T c3d = U11 * U10 * Zd1Q - Ud1 * Ud0 * Z11Q;
    T c2d = Ud0 * Z11S - U10 * Zd1S;

    T l3ZZASZ12 = l3_num * ZZAS * Z12;
    T l2Z12 = l2_num * Z12;

    T Zd2 = Zd1 * Z11Q * Z12 * ZALL * d;
    T Vd1 = l3ZZASZ12 * c3 + l2Z12 * c2 - V11 * c1;
    T Vd0 = l3ZZASZ12 * c3d + l2Z12 * c2d - V10 * c1;

    WeightedProjectiveMumford<T> ret(f, h, Ud1, Ud0, Vd1, Vd0, Zd1, Zd2);
    return ret;
}

template <class T>
WeightedProjectiveMumford<T> WeightedProjectiveMumford<T>::inv(){
    Polynomial<T> f = this->f;
    Polynomial<T> h = this->h;
    T U1 = this->U1;
    T U0 = this->U0;
    T V1 = this->V1;
    T V0 = this->V0;
    T Z1 = this->Z1;
    T Z2 = this->Z2;

    T h2 = this->h.coeff[2];
    T h1 = this->h.coeff[1];
    T h0 = this->h.coeff[0];

    T V1d = h2 * U1 * Z1 * Z2 - V1 - h1 * Z1 * Z1 * Z1 * Z2;
    T V0d = h2 * U0 * Z1 * Z2 - V0 - h0 * Z1 * Z1 * Z1 * Z2;
    WeightedProjectiveMumford<T> inv(f, h, U1, U0, V1d, V0d, Z1, Z2);
    return inv;
}

template <class T>
WeightedProjectiveMumford<T> WeightedProjectiveMumford<T>::zero(){
    Polynomial<T> f = this->f;
    Polynomial<T> h = this->h;

    WeightedProjectiveMumford<T> zero(f, h);
    return zero;
}

template <class T>
void WeightedProjectiveMumford<T>::print(){
    std::cout << "[";
    std::cout << this->U1 << " " << this->U0;
    std::cout << ", ";
    std::cout << this->V1 << " " << this->V0;
    std::cout << ", ";
    std::cout << this->Z1 << " " << this->Z2;
    std::cout << "]" << std::endl;

    Polynomial<T> u(2); Polynomial<T> v(1);
    T Z1S = this->Z1 * this->Z1;
    u.coeff[2] = T::ONE(); u.coeff[1] = this->U1 / Z1S; u.coeff[0] = this->U0 / Z1S;
    v.coeff[1] = this->V1 / (Z1S * this->Z1 * this->Z2); v.coeff[0] = this->V0 / (Z1S * this->Z1 * this->Z2);

    std::cout << "[" << u << ", " << v << "]" << std::endl;
    return;
}

template class WeightedProjectiveMumford<Number>;
