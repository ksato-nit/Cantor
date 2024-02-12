#include "mumford_projective.hpp"

template <class T>
ProjectiveMumford<T>::ProjectiveMumford(){
    this->f = Polynomial<T>();
    this->h = Polynomial<T>();
    this->U1 = T::ZERO();
    this->U0 = T::ONE();
    this->V1 = T::ZERO();
    this->V0 = T::ZERO();
    this->Z = T::ONE();
}

template <class T>
ProjectiveMumford<T>::ProjectiveMumford(Polynomial<T> f, Polynomial<T> h){
    this->f = f;
    this->h = h;
    this->U1 = T::ZERO();
    this->U0 = T::ONE();
    this->V1 = T::ZERO();
    this->V0 = T::ZERO();
    this->Z = T::ONE();
}

template <class T>
ProjectiveMumford<T>::ProjectiveMumford(Polynomial<T> f, Polynomial<T> h, T U1, T U0, T V1, T V0, T Z, T W1, T W0){
    this->f = f;
    this->h = h;
    this->U1 = U1;
    this->U0 = U0;
    this->V1 = V1;
    this->V0 = V0;
    this->Z = Z;
    this->W1 = W1;
    this->W0 = W0;
}

template <class T>
ProjectiveMumford<T>::ProjectiveMumford(Polynomial<T> f, Polynomial<T> h, T U1, T U0, T V1, T V0, T Z){
    this->f = f;
    this->h = h;
    this->U1 = U1;
    this->U0 = U0;
    this->V1 = V1;
    this->V0 = V0;
    this->Z = Z;
    this->W1 = U1 * U1;
    this->W0 = U1 * U0;
}

template <class T>
ProjectiveMumford<T>::ProjectiveMumford(Polynomial<T> f, Polynomial<T> h, T U1, T U0, T V1, T V0){
    this->f = f;
    this->h = h;
    this->U1 = U1;
    this->U0 = U0;
    this->V1 = V1;
    this->V0 = V0;
    this->Z = T::ONE();
    this->W1 = U1 * U1;
    this->W0 = U1 * U0;
}

template <class T>
ProjectiveMumford<T> ProjectiveMumford<T>::operator + (const ProjectiveMumford<T>& m) const{
    // deg u1 = deg u2 = 2
    ProjectiveMumford<T> ret = this->CostelloAdd(m);
    return ret;
}

template <class T>
ProjectiveMumford<T> ProjectiveMumford<T>::operator * (const mpz_class& k_) const{
    Polynomial<T> f = this->f;
    Polynomial<T> h = this->h;
    mpz_class k = k_;

    // double-and-add method によりスカラー倍を計算する．

    // 連除法で 2 進数に変換．
    mpz_class two = 2;
    int count = 0;
    std::vector<int> bits;
    while(true){
        mpz_class rem = k % 2;
        bits.push_back((int) rem.get_si());

        k = k / 2;
        ++count;

        if(k < 1){
            break;
        }
    }

    ProjectiveMumford<T> D = ProjectiveMumford<T>::zero(f, h);
    ProjectiveMumford<T> now = *this;
    for(int i = 0; i < count; ++i){
        if(bits[i] == 1){
            if(D.isZero()){
                D = now;
            }else if(now.isZero()){
                // D = D;
            }else{
                D = D.LangeAdd(now);
            }
            //D.print();
        }
        now = now.LangeDoubling();
    }
    return D;
}

template <class T>
ProjectiveMumford<T> ProjectiveMumford<T>::CostelloAdd(const ProjectiveMumford<T>& m) const{
    //std::cout << "Projective Costello Addition." << std::endl;

    // Total: 76M, 9S

    T U11 = this->U1;
    T U10 = this->U0;
    T U21 = m.U1;
    T U20 = m.U0;

    T V11 = this->V1;
    T V10 = this->V0;
    T V21 = m.V1;
    T V20 = m.V0;

    T Z1 = this->Z;
    T Z2 = m.Z;

    T f6 = this->f.coeff[6];
    T f5 = this->f.coeff[5];
    T f4 = this->f.coeff[4];

    // 10M, 3S
    T ZZ = Z1 * Z2;
    T U11Z2 = U11 * Z2;
    T U21Z1 = U21 * Z1;
    T U10Z2 = U10 * Z2;
    T U20Z1 = U20 * Z1;
    T V11Z2 = V11 * Z2;
    T V21Z1 = V21 * Z1;
    T V10Z2 = V10 * Z2;
    T V20Z1 = V20 * Z1;

    T U11Z2S = U11Z2 * U11Z2;
    T U21Z1S = U21Z1 * U21Z1;

    T T11 = U11 * U11;

    T Z1S = Z1 * Z1;
    T ZZ2 = ZZ * ZZ;

    // 3M

    // (Z1Z2)^2 がかかっている．
    T M1 = U11Z2S - U21Z1S + ZZ * (U20Z1 - U10Z2);
    T M2 = U21Z1 * U20Z1 - U11Z2 * U10Z2;
    // Z1Z2 がかかっている．
    T M3 = U11Z2 - U21Z1;
    T M4 = U20Z1 - U10Z2;

    T z1 = V10Z2 - V20Z1;
    T z2 = V11Z2 - V21Z1;

    // 4M

    // (Z1Z2)^3 がかかっている．
    T t1 = (M2 - z1) * (z2 - M1);
    T t2 = (-z1 - M2) * (z2 + M1);
    // (Z1Z2)^2 がかかっている．
    T t3 = (-z1 + M4) * (z2 - M3);
    T t4 = (-z1 - M4) * (z2 + M3);

    // 1M
    // (Z1Z2)^3 がかかっている．
    T l2_num = t1 - t2;
    // (Z1Z2)^2 がかかっている．
    T l3_num = t3 - t4;
    // (Z1Z2)^3 がかかっている．
    T d = -(M2 - M4) * (M1 + M3);
    d = d + d + (t3 + t4) - t1 - t2;

    // std::cout << l2_num << " " << l3_num << " " << d << " " << ZZ << std::endl;

    // 2M, 2S
    T A = d * d;
    T B = ZZ2 * l3_num * l3_num - f6 * A;

    // 6M
    // (l3_num^2 ZZ^2 - f6 d^2) ZZ^2 がかかっている．
    T l2l3ZZ = l2_num * l3_num * ZZ;
    T U1d = - B * (U11Z2 + U21Z1) - (f5 * A - l2l3ZZ - l2l3ZZ) * ZZ;
    T Zd = B * ZZ;
    
    // 21M, 1S
    T U0d = l3_num * ZZ2 * (U10 * Z1 - T11) + ZZ * Z1 * (l2_num * U11 + d * V11);
    U0d = (U0d + U0d) * l3_num + Z1S * (l2_num * l2_num - A * f4);
    U0d = Zd * Z2 * U0d;
    U0d = U0d - B * Z1 * ((U10Z2 + U20Z1 + U11 * U21) * Zd + (U11Z2 + U21Z1) * U1d);
    T ZdM = Z2 * Z1S * B;
    Zd = Zd * ZdM;
    U1d = U1d * ZdM;
    T ZdS = Zd * Zd;

    // 23M, 1S
    T dZ1ZdS = d * Z1 * ZdS;
    T ZZL3 = ZZ * l3_num;
    T Z1ZdL2 = Z1 * Zd * l2_num;
    T V1d = -Z1ZdL2 * (U11 * Zd - U1d * Z1) - ZZL3 * ( (U1d * U1d - U0d * Zd) * Z1S - (T11 - U10 * Z1) * ZdS) - V11 * dZ1ZdS;
    T V0d = Z1ZdL2 * (U0d * Z1 - U10 * Zd) + ZZL3 * (U11 * U10 * ZdS - U1d * U0d * Z1S) - V10 * dZ1ZdS;

    // 5M
    ZdM = Zd * Z1S * d;
    Zd = Zd * ZdM;

    U1d = U1d * ZdM;
    U0d = U0d * ZdM;

    ProjectiveMumford<T> ret(f, h, U1d, U0d, V1d, V0d, Zd);
    return ret;
}

template <class T>
ProjectiveMumford<T> ProjectiveMumford<T>::LangeAdd(const ProjectiveMumford<T>& m) const{
    //std::cout << "Projective Lange Addition." << std::endl;

    T U11 = this->U1;
    T U10 = this->U0;
    T U21 = m.U1;
    T U20 = m.U0;

    T V11 = this->V1;
    T V10 = this->V0;
    T V21 = m.V1;
    T V20 = m.V0;

    T Z1 = this->Z;
    T Z2 = m.Z;

    T f6 = this->f.coeff[6];
    T f5 = this->f.coeff[5];
    T f4 = this->f.coeff[4];

    // 64M, 6S

    // 1. 終結式を計算．
    // 8M, 2S
    T z1 = U11 * Z2 - U21 * Z1;
    T z2 = U20 * Z1 - U10 * Z2;
    T z3 = U11 * z1 + z2 * Z1;
    T r = z2 * z3 + z1 * z1 * U10; // Z1^3 Z2^2 がかかっている．
    T rs = r * r;

    // 2. almost inverse を計算．
    // 6M
    T inv1 = z1;
    T inv0 = z3;

    T w0 = V10 * Z2 - V20 * Z1;
    T w1 = V11 * Z2 - V21 * Z1;
    T w2 = inv0 * w0;
    T w3 = inv1 * w1;

    // 3. s を計算．
    // 4M
    T s1 = (inv0 + Z1 * inv1) * (w0 + w1) - w2 - w3 * (Z1 + U11);
    T s0 = w2 - U10 * w3;

    // 4. l を計算．全体に (Z1 Z2)^3 がかかっている．
    // 5M
    T l3 = s1 * Z2;
    T l2 = s1 * U21;
    T l0 = s0 * U20;
    T l1 = (s1 + s0) * (U21 + U20) - l2 - l0; //s1 * U20 + s0 * U21;
    l2 = l2 + s0 * Z2;

    // 5. U' を計算．全体に (Z1 Z2)^6 がかかっている．
    // 18M, 1S
    T f5Z2 = f5 * Z2;
    T f6U21 = f6 * U21;
    T rV21 = r * V21;

    T t4 = (s1 * l3 - rs * Z2 * f6) * Z2;
    T t3 = ((l2 * s1 + l3 * s0) - rs * (f5Z2 - f6U21)) * Z2;
    T t2 = Z2 * (s0 * l2 + s1 * (l1 + rV21 * 2)) - rs * ( (f4 * Z2 - f6 * U20) * Z2 - (f5Z2 - f6U21) * U21 );
 
    // 8M, 2S
    T t4U11 = t4 * U11;
    T t3Z1 = t3 * Z1;
    T Ud2 = t4;
    Ud2 = Ud2 * Z1 * Z1;
    T Ud1 = t3Z1 - t4U11;
    Ud1 = Ud1 * Z1;
    T Ud0 = (t2 * Z1 - t4 * U10) * Z1 - (t3Z1 - t4U11) * U11;
    T Zd = Ud2;
    T ZdS = Zd * Zd;

    // 6. V' を計算．
    // 10M, 1S
    T Vd0 = (-l0 - V20 * r) * ZdS - Ud0 * (Ud1 * l3 - l2 * Zd);
    T Vd1 = (-l1 - rV21) * ZdS - l3 * (Ud1 * Ud1 - Ud0 * Zd) + Zd * Ud1 * l2;

    // 7. Z' を調整．
    // 5M
    T M = Zd * Z2 * r;
    Ud1 = Ud1 * M;
    Ud0 = Ud0 * M;
    Zd = Zd * M;

    ProjectiveMumford<T> ret(f, h, Ud1, Ud0, Vd1, Vd0, Zd);
    return ret;
}

template <class T>
ProjectiveMumford<T> ProjectiveMumford<T>::LangeDoubling() const{
    //std::cout << "Projective Lange Doubling." << std::endl;

    T U1 = this->U1;
    T U0 = this->U0;

    T V1 = this->V1;
    T V0 = this->V0;

    T Z = this->Z;

    T f6 = this->f.coeff[6];
    T f5 = this->f.coeff[5];
    T f4 = this->f.coeff[4];
    T f3 = this->f.coeff[3];
    T f2 = this->f.coeff[2];

    // 59M, 9S

    // 1. precomputation.
    // 3M, 2S

    T U0Z = U0 * Z;
    T f5Z = f5 * Z;
    T Z2 = Z * Z;
    T Z3 = Z2 * Z;
    T Z4 = Z2 * Z2;

    // 2. v~ と u の終結式を計算．
    // 4M, 2S
    T V1t = V1 * 2;
    T V0t = V0 * 2;

    T W0 = V1 * V1;
    T U1s = U1 * U1;
    T W1 = U1s;
    T W2 = W0 * 4;
    T W3 = U1 * V1t;
    // Z^3 がかかっている．
    T V0tZ = V0t * Z;
    T r = U0 * W2 + V0t * (V0tZ - W3);

    // 3. almost inverse を計算．
    // Z がかかっている．
    T inv1d = -V1t;
    // Z^2 がかかっている．
    T inv0d = V0tZ - W3;

    // 4. k を計算．
    // 15M
    T k4 = f6;
    T k4U1 = k4 * U1;
    T k4U0Z = k4 * U0Z;
    T k3 = f5Z - k4U1;
    T k3U0Z = k3 * U0Z;
    T f4Z2 = f4 * Z2;
    T k2 = f4Z2 - k4U0Z - k3 * U1;
    T k1 = f3 * Z3 - k3U0Z - k2 * U1;
    T k0 = f2 * Z4 - W0 * Z2 - k2 * U0Z - k1 * U1;
    // Z^3 がかかっている．
    T k1d = k1 + W1 * (k3 - k4U1) - k3U0Z + U1 * (k4U0Z * 2 - k2);
    // Z^4 がかかっている．
    T k0d = k0 + U0Z * (U1 * (k3 - k4U1) + k4U0Z - k2);

    // 5. s を計算．
    // 5M
    W0 = k0d * inv0d;
    W1 = k1d * inv1d;
    // Z^5 がかかっている．
    T s1d = (inv0d + inv1d) * (k0d + k1d) - W0 - W1 * (T::ONE() + U1);
    // Z^6 がかかっている．
    T s0d = W0 - W1 * U0Z;

    // 6. U' を計算．
    // 12M, 4S
    T rs = r * r;
    T f5Zk = f5Z - k4U1 * 2;
    T Z3r = Z3 * r;
    T rsZ4 = rs * Z4;
    // Z^10 がかかっている．
    T Ud2 = s1d * s1d - rsZ4 * f6;
    // Z^11 がかかっている．
    T Ud1 = s1d * s0d * 2 - rsZ4 * f5Zk;
    // Z^12 がかかっている．
    T V1Z3r = V1 * Z3r;
    T Ud0 = s0d * s0d + s1d * V1Z3r * 2 - rsZ4 * (f4Z2 - (U0Z * 2 + U1s) * f6 - U1 * f5Zk * 2);
    Ud1 = Ud1 * Z;
    T Zd = Ud2 * Z2;
    T Zd2 = Zd * Zd;

    // 7. V' を計算．
    // 15M, 1S
    // r Z^7 がかかっている．
    T l3 = s1d * Z2;
    T l2 = (s1d * U1 + s0d) * Z;
    T l1 = (s1d * U0Z + s0d * U1);
    T l0 = s0d * U0;
    T l2Zd = l2 * Zd;
    T Vd1 = -l3 * (Ud1 * Ud1 - Ud0 * Zd) + Ud1 * l2Zd - (l1 + V1Z3r) * Zd2;
    T Vd0 = -(l3 * Ud1 - l2Zd) * Ud0 - (l0 + V0 * Z3r) * Zd2;

    // 8. 調整
    // 5M
    T ZM = Z4 * Zd * r;
    Ud1 = Ud1 * ZM;
    Ud0 = Ud0 * ZM;
    Zd = Zd * ZM;


    ProjectiveMumford<T> ret(f, h, Ud1, Ud0, Vd1, Vd0, Zd);
    return ret;
}

template <class T>
ProjectiveMumford<T> ProjectiveMumford<T>::inv(){
    Polynomial<T> f = this->f;
    Polynomial<T> h = this->h;
    T U1 = this->U1;
    T U0 = this->U0;
    T V1 = this->V1;
    T V0 = this->V0;
    T Z = this->Z;

    T h2 = this->h.coeff[2];
    T h1 = this->h.coeff[1];
    T h0 = this->h.coeff[0];
    // - (h + v) % u を計算．
    T V1d = U1 * h2 - (V1 + h1) * Z;
    T V0d = U0 * h2 - (V0 + h0) * Z;
    ProjectiveMumford<T> inv(f, h, U1, U0, V1d, V0d, Z);
    return inv;
}

template <class T>
ProjectiveMumford<T> ProjectiveMumford<T>::zero() const{
    Polynomial<T> f = this->f;
    Polynomial<T> h = this->h;

    ProjectiveMumford<T> zero(f, h);
    return zero;
}

template <class T>
ProjectiveMumford<T> ProjectiveMumford<T>::zero(const Polynomial<T>& f, const Polynomial<T>& h){
    ProjectiveMumford<T> zero(f, h);
    return zero;
}

template <class T>
void ProjectiveMumford<T>::print() const{
    std::cout << "[";
    std::cout << this->U1 << " " << this->U0;
    std::cout << ", ";
    std::cout << this->V1 << " " << this->V0;
    std::cout << ", ";
    std::cout << this->Z;
    std::cout << "]" << std::endl;

    Polynomial<T> u(2); Polynomial<T> v(1);
    u.coeff[2] = T::ONE(); u.coeff[1] = this->U1 / this->Z; u.coeff[0] = this->U0 / this->Z;
    v.coeff[1] = this->V1 / this->Z; v.coeff[0] = this->V0 / this->Z;

    std::cout << "[" << u << ", " << v << "]" << std::endl;
    
    return;
}

template <class T>
bool ProjectiveMumford<T>::isZero() const{
    T U1 = this->U1;
    T U0 = this->U0;
    T V1 = this->V1;
    T V0 = this->V0;
    T Z = this->Z;

    if(V1.isZero() && V0.isZero()){
        if(U1.isZero() && U0 == Z){
            return true;
        }
    }
    return false;
}

template class ProjectiveMumford<Number>;
template class ProjectiveMumford<ExtendedNumber>;
