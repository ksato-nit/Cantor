#include "mumford.hpp"

template <class T>
Mumford<T>::Mumford(){
    this->f = Polynomial<T>();
    this->h = Polynomial<T>();
    this->u2 = T::ONE();
    this->u1 = T::ZERO();
    this->u0 = T::ONE();
    this->v1 = T::ZERO();
    this->v0 = T::ZERO();
}

template <class T>
Mumford<T>::Mumford(Polynomial<T> f, Polynomial<T> h){
    this->f = f;
    this->h = h;
    this->u2 = T::ONE();
    this->u1 = T::ZERO();
    this->u0 = T::ONE();
    this->v1 = T::ZERO();
    this->v0 = T::ZERO();
}

template <class T>
Mumford<T>::Mumford(Polynomial<T> f, Polynomial<T> h, T u1, T u0, T v1, T v0){
    this->f = f;
    this->h = h;
    this->u2 = T::ONE();
    this->u1 = u1;
    this->u0 = u0;
    this->v1 = v1;
    this->v0 = v0;
}

// todo: 書き換え
template <class T>
Mumford<T> Mumford<T>::CostelloScalarMultiple (const mpz_class& k_) const{
    Polynomial<T> f = this->f;
    Polynomial<T> h = this->h;
    mpz_class k = k_;

    // operator * の，Costello による計算．
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

    Mumford<T> D = Mumford<T>::zero(f, h);
    Mumford<T> now = *this;
    for(int i = 0; i < count; ++i){
        if(bits[i] == 1){
            if(D.isZero()){
                D = now;
            }else if(now.isZero()){
                // D = D;
            }else{
                D = D.CostelloAdd(now);
            }
            //D.print();
        }
        now = now.CostelloDoubling();
    }
    return D;
}

template <class T>
Mumford<T> Mumford<T>::operator * (const mpz_class& k_) const{
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

    Mumford<T> D = Mumford<T>::zero(f, h);
    Mumford<T> now = *this;
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
Mumford<T> Mumford<T>::operator + (const Mumford<T>& m) const{
    if(this->f.deg == 6){
        return this->LangeAdd(m);
    }else{
        return this->CantorAdd(m);
    }
}

template <class T>
Mumford<T> Mumford<T>::CostelloAdd(const Mumford<T>& m) const{
    //std::cerr << "Costello Addition." << std::endl;
    T u11 = this->u1;
    T u10 = this->u0;
    T u21 = m.u1;
    T u20 = m.u0;

    T v11 = this->v1;
    T v10 = this->v0;
    T v21 = m.v1;
    T v20 = m.v0;

    T f6;
    if(this->f.deg == 6){
        f6 = this->f.coeff[6];
    }else{
        f6 = T::ZERO();
    }
    T f5 = this->f.coeff[5];
    T f4 = this->f.coeff[4];

    // TODO: U1, U0 を座標に含めて保持するようにする．
    T U11 = u11 * u11;
    T U21 = u21 * u21;
    T U10 = u11 * u10;
    T U20 = u21 * u20;

    T u1S = u11 + u21;
    T v0D = v10 - v20;
    T v1D = v11 - v21;

    T M1 = U11 - U21 - u10 + u20;
    T M2 = U20 - U10;
    T M3 = u11 - u21;
    T M4 = u20 - u10;

    T t1 = (M2 - v0D) * (v1D - M1);
    T t2 = -(v0D + M2) * (v1D + M1);
    T t3 = (M4 - v0D) * (v1D - M3);
    T t4 = -(v0D + M4) * (v1D + M3);

    T l2_num = t1 - t2;
    T l3_num = t3 - t4;

    T d = (M4 - M2) * (M1 + M3);
    d = d + d;
    d = d + t3 + t4 - t1 - t2;

    T A = d * d;
    T B = l3_num * l3_num - f6 * A;
    T C = (d * B).inv();
    T d_inv = B * C;
    T d_shifted_inv = d * A * C;

    T l2 = l2_num * d_inv;
    T l3 = l3_num * d_inv;

    T l0 = v10 + l2 * u10 - l3 * U10;
    T l1 = v11 + l2 * u11 - l3 * (U11 - u10);

    T u1dd = -(u1S + (f5 - l2 * l3 - l2 * l3) * d_shifted_inv);

    T u0dd = l3 * (l3 * (u10 - U11) + l2 * u11 + v11);
    u0dd = u0dd + u0dd;
    u0dd = u0dd + l2 * l2 - f4;
    u0dd = u0dd * d_shifted_inv;
    u0dd = u0dd - u11 * u21 - u10 - u20 - u1S * u1dd;

    T U1dd = u1dd * u1dd;
    T U0dd = u1dd * u0dd;

    T v1dd = l3 * (u0dd - U1dd + U11 - u10) + l2 * (u1dd - u11) - v11;
    T v0dd = l3 * (U10 - U0dd) + l2 * (u0dd - u10) - v10;

    return Mumford<T>(f, h, u1dd, u0dd, v1dd, v0dd);
}

template <class T>
Mumford<T> Mumford<T>::HarleyAdd(const Mumford<T>& m) const{
    std::cerr << "Harley Addition." << std::endl;
    T u11 = this->u1;
    T u10 = this->u0;
    T u21 = m.u1;
    T u20 = m.u0;

    T v11 = this->v1;
    T v10 = this->v0;
    T v21 = m.v1;
    T v20 = m.v0;

    T h0 = this->h.coeff[0];
    T h1 = this->h.coeff[1];
    T h2 = this->h.coeff[2];

    T f4 = this->f.coeff[4];

    // 1. u1, u2 の終結式を計算．
    T z1 = u11 - u21;
    T z2 = u20 - u10;
    T z3 = u11 * z1 + z2;
    T r = z2 * z3 + z1 * z1 * u10;

    // 2. u2 の almost inverse (mod u1) を計算．
    T inv1 = z1;
    T inv0 = z3;

    // 3. s' を計算．
    T w0 = v10 - v20;
    T w1 = v11 - v21;
    T w2 = inv0 * w0;
    T w3 = inv1 * w1;
    T s1d = (inv0 + inv1) * (w0 + w1) - w2 - w3 * (T::ONE() + u11);
    T s0d = w2 - u10 * w3;

    if(!s1d.isZero()){
        // 4. s'' を計算．
        T w1 = (r * s1d).inv();
        T w2 = r * w1;
        T w3 = s1d * s1d * w1;
        T w4 = r * w2;
        T w5 = w4 * w4;
        T s0dd = s0d * w2;

        // 5. l' を計算．
        T l2d = u21 + s0dd;
        T l1d = u21 * s0dd + u20;
        T l0d = u20 * s0dd;


        // 6. u' を計算．
        T u0d = (s0dd - u11) * (s0dd - z1 + h2 * w4) - u10 + l1d + (h1 + v21 * 2) * w4 + (u21 * 2 + z1 - f4) * w5;
        T u1d = s0dd * 2 - z1 + h2 * w4 - w5;

        // 7. v' を計算．
        w1 = l2d - u1d;
        w2 = u1d * w1 + u0d - l1d;
        T v1d = w2 * w3 - v21 - h1 + h2 * u1d;
        w2 = u0d * w1 - l0d;
        T v0d = w2 * w3 - v20 - h0 + h2 * u0d;

        Mumford ret(f, h, u1d, u0d, v1d, v0d);
        return ret;
    }else{
        //std::cerr << "Special case." << std::endl;
        // サブルーチン

        // 4'. s を計算．
        T inv = r.inv();
        T s0 = s0d * inv;

        // 5'. u' を計算．
        T u0d = f4 - u21 - u11 - s0*s0 - s0 * h2;

        // 6'. v' を計算．
        T w1 = s0 * (u21 + u0d) + h1 + v21 + h2 * u0d;
        T w2 = s0 + v20 + h0;
        T v0d = u0d * w1 - w2;

        Mumford<T> ret(f, h, T::ONE(), u0d, T::ZERO(), v0d);
        return ret;
    }
}

template <class T>
Mumford<T> Mumford<T>::HarleyAddDegenerated(const Mumford<T>& m) const{
    T u11 = this->u1;
    T u10 = this->u0;
    T u21 = m.u1;
    T u20 = m.u0;

    T v11 = this->v1;
    T v10 = this->v0;
    T v21 = m.v1;
    T v20 = m.v0;

    T h0 = this->h.coeff[0];
    T h1 = this->h.coeff[1];
    T h2 = this->h.coeff[2];

    T f3 = this->f.coeff[3];
    T f4 = this->f.coeff[4];

    // 1. r を計算．
    T r = u20 - (u21 - u10) * u10;

    // 2. u2 の almost inverse (mod u1) を計算．
    T inv = r.inv();

    // 3. s を計算．
    T s0 = inv * (v10 - v20 + v21 * u10);

    // 4. l を計算．
    T l1 = s0 * u21;
    T l0 = s0 * u20;

    // 5. k を計算．
    T k2 = f4 - u21;
    T k1 = f3 - (f4 - u21) * u21 - v21 * h2 - u20;

    // 6. u' を計算．
    T u1d = k2 - s0 * s0 - s0 * h2 - u10;
    T u0d = k1 - s0 * (l1 + h1 + v21 * 2) - u10 * u1d;

    // 7. v' を計算．
    T v1d = (h2 + s0) * u1d - (h1 + l1 + v21);
    T v0d = (h2 + s0) * u0d - (h0 + l0 + v20);

    Mumford<T> ret(f, h, u1d, u0d, v1d, v0d);
    return ret;
}

template <class T>
Mumford<T> Mumford<T>::LangeAdd(const Mumford<T>& m) const{
    //std::cerr << "Lange Addition." << std::endl;
    T u11 = this->u1;
    T u10 = this->u0;
    T u21 = m.u1;
    T u20 = m.u0;

    T v11 = this->v1;
    T v10 = this->v0;
    T v21 = m.v1;
    T v20 = m.v0;

    T f6 = this->f.coeff[6];
    T f5 = this->f.coeff[5];
    T f4 = this->f.coeff[4];

    // 1. u1, u2 の終結式を計算．
    T z1 = u11 - u21;
    T z2 = u20 - u10;
    T z3 = u11 * z1 + z2;
    T r = z2 * z3 + z1 * z1 * u10;

    // 2. u2 の almost inverse (mod u1) を計算．
    T inv1 = z1;
    T inv0 = z3;

    // 3. s' を計算．
    T w0 = v10 - v20;
    T w1 = v11 - v21;
    T w2 = inv0 * w0;
    T w3 = inv1 * w1;
    T s1d = (inv0 + inv1) * (w0 + w1) - w2 - w3 * (T::ONE() + u11);
    T s0d = w2 - u10 * w3;

    if(s1d.isZero()){
        //std::cerr << "Special case." << std::endl;
        // todo: ここの場合分けを厳密に書く．
        Mumford ret(f, h);
        return ret;
    }

    // 4. l' を計算．
    T l3d = s1d;
    T l2d = u21 * s1d;
    T l0d = u20 * s0d;
    T l1d = (s1d + s0d) * (u21 + u20) - l2d - l0d; //u21 * s0d + u20 * s1d; 
    l2d = l2d + s0d;

    // 5. u' を計算．
    T k4 = f6;
    T k3 = f5 - f6 * u21;
    T k2 = f4 - f6 * u20 - (f5 - f6 * u21) * u21;

    T t4 = s1d * l3d - k4 * r * r;
    T t3 = s1d * l2d + s0d * l3d - k3 * r * r;
    T t2 = r * v21;
    t2 = t2 + t2;
    t2 = s1d * (t2 + l1d) + s0d * l2d - k2 * r * r;

    T u0d = t2 - t4 * u10 - (t3 - t4 * u11) * u11;
    T u1d = t3 - t4 * u11;
    T u2d = t4;

    // 6. u' をモニックにする．
    w1 = (u2d * r).inv();
    w2 = w1 * r;
    w3 = w1 * u2d;
    u1d = u1d * w2;
    u0d = u0d * w2;
    // ud の計算まで正しい．

    // 7. v' を計算．
    T v1d = (-l1d + (u0d - u1d * u1d) * l3d + u1d * l2d) * w3 - v21;
    T v0d = (-l0d  - u1d * u0d * l3d + u0d * l2d) * w3 - v20;

    Mumford<T> ret(f, h, u1d, u0d, v1d, v0d);

    return ret;
}

template <class T>
Mumford<T> Mumford<T>::LangeDoubling() const{
    //std::cerr << "Lange Doubling." << std::endl;
    T u1 = this->u1;
    T u0 = this->u0;
    T v1 = this->v1;
    T v0 = this->v0;

    T f2 = this->f.coeff[2];
    T f3 = this->f.coeff[3];
    T f4 = this->f.coeff[4];
    T f5 = this->f.coeff[5];
    T f6 = this->f.coeff[6];

    // 44M, 6S, I

    // 1. v~ (=2v) と u の終結式を計算．
    // 3M, 2S
    T v1t = v1 + v1;
    T v0t = v0 + v0;

    T w0 = v1 * v1;
    T u1s = u1 * u1;
    T w1 = u1s;
    T w2 = w0 + w0;
    w2 = w2 + w2;
    T w3 = u1 * v1t;
    T r = u0 * w2 + v0t * (v0t - w3);

    // 2. r の almost inverse を計算．
    T inv1d = -v1t;
    T inv0d = v0t - w3;

    // 3. k を計算．
    // 11M
    T k4 = f6;
    T k4u0 = k4 * u0;
    T k4u1 = k4 * u1;
    T k3 = f5 - k4u1;
    T k3u0 = k3 * u0;
    T k2 = f4 - k4u0 - k3 * u1;
    T k1 = f3 - k3u0 - k2 * u1;
    T k0 = f2 - w0 - k2 * u0 - k1 * u1;
    T u1kd = u1 * (k3 - k4u1);
    T k1d = k1 + w1 * (k3 - k4u1) - k3u0 + u1 * (k4u0 + k4u0 - k2);
    T k0d = k0 + u0 * (u1kd + k4u0 - k2);

    // 4. s' を計算．
    // 5M
    w0 = k0d * inv0d;
    w1 = k1d * inv1d;
    T s1d = (inv0d + inv1d) * (k0d + k1d) - w0 - w1 * (T::ONE() + u1);
    T s0d = w0 - u0 * w1;

    if(s1d.isZero()){
        std::cerr << "Special case." << std::endl;
        // サブルーチン
        Mumford ret(f, h);
        return ret;
    }

    // 5. u' を計算．
    // 8M, 3S
    T rs = r * r;
    T u2d = s1d * s1d - rs * f6;
    T u1d = s1d * s0d;
    u1d = u1d + u1d;
    u1d = u1d - rs * (f5 - k4u1 - k4u1);
    T u0d = s0d * s0d + v1 * s1d * r * 2 - rs * (f4 - (u0 + u0 + u1s) * f6 - u1 * (f5 - k4u1 - k4u1) * 2);

    // 6. r と u2d の逆元を計算．
    // 3M, I
    w0 = (r * u2d).inv();
    w1 = w0 * r;
    w2 = w0 * u2d;

    // 7. u' を計算．
    // 2M
    u1d = u1d * w1;
    u0d = u0d * w1;

    // 8. v' を計算．
    // 12M, 1S
    T l3 = s1d * w2;
    T l2 = (s1d * u1 + s0d) * w2;
    T l1 = (s1d * u0 + s0d * u1) * w2;
    T l0 = s0d * u0 * w2;
    T v1d = l3 * (u1d * u1d - u0d) - l2 * u1d + l1 + v1;
    T v0d = (l3 * u1d - l2) * u0d + l0 + v0;

    Mumford<T> ret(f, h, u1d, u0d, -v1d, -v0d);
    return ret;
}

template <class T>
Mumford<T> Mumford<T>::CostelloDoubling() const{
    //std::cerr << "Costello Doubling." << std::endl;
    T u1 = this->u1;
    T u0 = this->u0;
    T v1 = this->v1;
    T v0 = this->v0;

    T U1 = u1 * u1;
    T U0 = u1 * u0;

    T f2 = this->f.coeff[2];
    T f3 = this->f.coeff[3];
    T f4 = this->f.coeff[4];
    T f5 = this->f.coeff[5];
    T f6 = this->f.coeff[6];

    // 32M, 6S, I

    T vv = v1 * v1;
    T va = (v1 + u1) * (v1 + u1) - vv - U1;

    T M1 = (v0 - va) * 2;
    T M2 = (U1 * 2 + u0) * v1 * 2;
    T M3 = -v1 * 2;
    T M4 = va + v0 * 2;

    T f6u0 = f6 * u0;
    T f6U1 = f6 * U1;
    T f5u0 = f5 * u0;
    T f5u1 = f5 * u1;

    T z11 = f5u1 * 2 - f6U1 * 3 - f4;
    z11 = z11 * U1;
    T z12 = f5u1 * 2 + f6u0 * 3 - f4 * 2;
    z12 = z12 * u0;

    T z1 = z11 + z12 - vv + f2;
    T z2 = (f6u0 * 6 - f6U1 * 4 + f5u1 * 3 - f4 * 2) * u1 - f5u0 * 2 + f3;

    T t11 = M2 - z1;
    T t12 = z2 - M1;
    T t21 = -(z1 + M2);
    T t22 = z2 + M1;
    T t31 = -z1 + M4;
    T t32 = z2 - M3;
    T t41 = -(z1 + M4);
    T t42 = z2 + M3;

    T t1 = t11 * t12;
    T t2 = t21 * t22;
    T t3 = t31 * t32;
    T t4 = t41 * t42;

    T l2_num = t1 - t2;
    T l3_num = t3 - t4;

    T d11 = M4 - M2;
    T d12 = M1 + M3;
    T d = d11 * d12 * 2 + t3 + t4 - t1 - t2;

    T A = d * d;
    T B1 = l3_num * l3_num;
    T B2 = f6 * A;
    T B = B1 - B2;
    T C = (d * B).inv();

    T d_inv = B * C;
    T d_shifted_inv = A * C * d;

    T l2 = l2_num * d_inv;
    T l3 = l3_num * d_inv;

    T u1d = (l2 * l3 * 2 - f5) * d_shifted_inv - u1 * 2;
    T u0d = (((u0 - U1) * l3 + l2 * u1 + v1) * l3 * 2 + l2 * l2 - f4) * d_shifted_inv;
    u0d = u0d - u0 * 2 - u1d * u1 * 2 - U1;
    T U1d = u1d * u1d;
    T U0d = u1d * u0d;

    T v1d = (u0d - U1d + U1 - u0) * l3 + (u1d - u1) * l2 - v1;
    T v0d = (U0 - U0d) * l3 + (u0d - u0) * l2 - v0;

    Mumford<T> ret(f, h, u1d, u0d, v1d, v0d);
    return ret;
}

template <class T>
Mumford<T> Mumford<T>::CantorAdd(const Mumford<T>& m) const{
    std::cerr << "Cantor Addition." << std::endl;

    Polynomial<T> h = this->h;

    Polynomial<T> u1 = Polynomial<T>(2);
    u1.coeff[2] = T::ONE();
    u1.coeff[1] = this->u1;
    u1.coeff[0] = this->u0;
    Polynomial<T> v1 = Polynomial<T>(1, this->v0, this->v1);
    Polynomial<T> u2 = Polynomial<T>(2);
    u2.coeff[2] = T::ONE();
    u2.coeff[1] = m.u1;
    u2.coeff[0] = m.u0;
    Polynomial<T> v2 = Polynomial<T>(1, m.v0, m.v1);

    auto tup1 = Polynomial<T>::extended_gcd(u1, u2);
    Polynomial<T> d1 = std::get<0>(tup1);
    Polynomial<T> e1 = std::get<1>(tup1);
    Polynomial<T> e2 = std::get<2>(tup1);

    auto tup2 = Polynomial<T>::extended_gcd(d1, v1 + v2 + h);
    Polynomial<T> d = std::get<0>(tup2);
    Polynomial<T> c1 = std::get<1>(tup2);
    Polynomial<T> c2 = std::get<2>(tup2);

    Polynomial<T> s1 = c1 * e1;
    Polynomial<T> s2 = c1 * e2;

    Polynomial<T> u = u1 * u2 / (d * d);
    Polynomial<T> v = ((s1 * u1 * v2 + s2 * u2 * v1 + c2 * (v1 * v2 + f)) / d) % u;

    while(u.deg > GENUS){
        Polynomial<T> ud = (f - (v * h) - (v * v)) / u;
        Polynomial<T> vr = (-(h + v)) % ud;

        u = ud;
        v = vr;
    }

    T lc = u.coeff[u.deg];
    u = u * lc.inv();

    T v1n = (v.deg == 1) ? v.coeff[1] : T::ZERO();

    Mumford<T> ret(f, h, u.coeff[1], u.coeff[0], v1n, v.coeff[0]);

    return ret;
}

template <class T>
Mumford<T> Mumford<T>::inv(){
    Polynomial<T> f = this->f;
    Polynomial<T> h = this->h;
    Polynomial<T> u(2);
    u.coeff[0] = u0;
    u.coeff[1] = u1;
    u.coeff[2] = u2;
    Polynomial<T> v(1, v0, v1);

    Polynomial<T> vd = (-(h + v)) % u;

    Mumford<T> inv(f, h, u1, u0, vd.coeff[1], vd.coeff[0]);
    return inv;
}

template <class T>
Mumford<T> Mumford<T>::zero() const{
    Polynomial<T> f = this->f;
    Polynomial<T> h = this->h;

    Mumford<T> zero(f, h);
    return zero;
}

template <class T>
Mumford<T> Mumford<T>::zero(const Polynomial<T>& f, const Polynomial<T>& h){
    Mumford<T> zero(f, h);
    return zero;
}

template <class T>
void Mumford<T>::print() const{
    Polynomial<T> u(2);
    u.coeff[0] = u0;
    u.coeff[1] = u1;
    u.coeff[2] = Number::ONE();
    Polynomial<T> v(1, v1, v0);
    std::cerr << "[" << u << ", " << v << "]" << std::endl;
    return;
}

template <class T>
bool Mumford<T>::isZero() const{
    if(v1.isZero() && v0.isZero()){
        if(u1.isZero() && u0 == Number::ONE()){
            return true;
        }
    }
    return false;
}

template class Mumford<Number>;
