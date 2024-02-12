#include "polynomial.hpp"

template <class T>
Polynomial<T>::Polynomial(){
    this->deg = 0;
    this->coeff.resize(1);
    this->coeff[0] = T::ZERO();
}

template <class T>
Polynomial<T>::Polynomial(int deg){
    this->deg = deg;
    this->coeff.resize(deg + 1);
}

template <class T>
Polynomial<T>::Polynomial(int deg, int* coeff){
    this->deg = deg;
    this->coeff.resize(deg + 1);
    for(int i = 0; i <= deg; ++i){
        T c(coeff[i]);
        this->coeff[i] = c;
    }
}

template <class T>
Polynomial<T>::Polynomial(std::vector<T> coeff){
    int deg = coeff.size() + 1;
    this->deg = deg;
    this->coeff.resize(deg + 1);
    for(int i = 0; i <= deg; ++i){
        this->coeff[i] = coeff[i];
    }
}

template <class T>
Polynomial<T>::Polynomial(int deg, int c0){
    this->deg = deg;
    this->coeff.resize(deg + 1);
    T c(c0);
    this->coeff[0] = c;
}

template <class T>
Polynomial<T>::Polynomial(int deg, T c0){
    this->deg = deg;
    this->coeff.resize(deg + 1);
    this->coeff[0] = c0;
}

template <class T>
Polynomial<T>::Polynomial(int deg, T c0, T c1){
    this->deg = deg;
    this->coeff.resize(deg + 1);
    this->coeff[0] = c0;
    this->coeff[1] = c1;
}

template <class T>
void Polynomial<T>::resize(int deg){
    this->deg = deg;
    this->coeff.resize(deg + 1);
}

template <class T>
void Polynomial<T>::resize(std::vector<T> coeff){
    this->deg = coeff.size() + 1;
    this->coeff = coeff;
}

// 次数が整合的であるようにする (高次の係数で 0 のものがあれば切り落としていく)．
template <class T>
void Polynomial<T>::normalize(){
    int deg = this->deg;
    for(auto i = this->coeff.rbegin(), e = this->coeff.rend(); i != e; ++i){
        if(i->isZero()){
            --deg;
        }else{
            break;
        }
    }

    if(deg <= 0){
        deg = 0;
    }

    this->deg = deg;
    this->coeff.resize(deg + 1);
}

template <class T>
bool Polynomial<T>::isZero(){
    for(auto n : this->coeff){
        if(!n.isZero()){
            return false;
        }
    }
    return true;
}

template <class T>
Polynomial<T> Polynomial<T>::operator + (const Polynomial<T>& q) const{
    // TODO : もう少しスマートな実装を考える．
    int deg = std::max(this->deg, q.deg);
    Polynomial<T> r(deg);
    for(int i = 0; i <= deg; ++i){
        T a = (this->deg >= i) ? this->coeff[i] : T::ZERO();
        T b = (q.deg >= i) ? q.coeff[i] : T::ZERO();
        r.coeff[i] = a + b;
    }
    //r.normalize();
    return r;
}

template <class T>
Polynomial<T> Polynomial<T>::operator - (const Polynomial<T>& q) const{
    int deg = std::max(this->deg, q.deg);
    Polynomial<T> r(deg);
    for(int i = 0; i <= deg; ++i){
        T a = (this->deg >= i) ? this->coeff[i] : T::ZERO();
        T b = (q.deg >= i) ? q.coeff[i] : T::ZERO();
        r.coeff[i] = a - b;
    }
    //r.normalize();
    return r;
}

template <class T>
Polynomial<T> Polynomial<T>::operator * (const Polynomial<T>& q) const{
    int deg = this->deg + q.deg;
    Polynomial<T> r(deg);
    for(int i = 0; i <= this->deg; ++i){
        for(int j = 0; j <= q.deg; ++j){
            T a = (this->deg >= i) ? this->coeff[i] : T::ZERO();
            T b = (q.deg >= j) ? q.coeff[j] : T::ZERO();
            r.coeff[i + j] = r.coeff[i + j] + a * b;
        }
    }
    //r.normalize();
    return r;
}

template <class T>
Polynomial<T> Polynomial<T>::operator * (const T& k) const{
    int deg = this->deg;
    Polynomial<T> r(deg);
    for(int i = 0; i <= this->deg; ++i){
        r.coeff[i] = this->coeff[i] * k;
    }
    return r;
}

template <class T>
Polynomial<T> Polynomial<T>::operator / (const Polynomial<T>& g) const{
    auto tup = Polynomial<T>::divide(*this, g);

    Polynomial<T> q = std::get<0>(tup);
    Polynomial<T> r = std::get<1>(tup);

    return q;
}

template <class T>
Polynomial<T> Polynomial<T>::operator % (const Polynomial<T>& g) const{
    if(this->deg < g.deg){
        return *this;
    }

    auto tup = Polynomial<T>::divide(*this, g);

    Polynomial<T> q = std::get<0>(tup);
    Polynomial<T> r = std::get<1>(tup);

    return r;
}

template <class T>
Polynomial<T> Polynomial<T>::operator + () const{
    return *this;
}

template <class T>
Polynomial<T> Polynomial<T>::operator - () const{
    Polynomial<T> r = (*this) * T::MINUS_ONE();
    return r;
}

template <class T>
bool Polynomial<T>::operator == (const Polynomial<T>& g) const{
    if(this->deg != g.deg){
        return false;
    }

    for(int i = 0; i <= this->deg; ++i){
        if(this->coeff[i] != g.coeff[i]){
            return false;
        }
    }

    return true;
}

template <class T>
bool Polynomial<T>::operator != (const Polynomial<T>& g) const{
    return !(*this == g);
}

// f = qg + r を満たす q, r を求める．
// deg g <= deg f を仮定．
template <class T>
std::tuple<Polynomial<T>, Polynomial<T>> Polynomial<T>::divide(Polynomial<T> f, Polynomial<T> g){
    int deg_f = f.deg;
    int deg_g = g.deg;
    Polynomial<T> q(deg_f - deg_g);
    Polynomial<T> r(deg_f);

    // TODO : インデックスをシフトしてもう少しきれいに書けそう．
    for(int j = deg_f - deg_g; j >= 0; --j){
        q.coeff[j] = f.coeff[j + deg_g] / g.coeff[deg_g];
        Polynomial<T> mono(j);
        mono.coeff[j] = q.coeff[j];
        r = f - (g * mono);
        f = r;
    }

    r.normalize();

    return std::forward_as_tuple(q, r);
}


// s = gcd(f, g), fu + gv = s を満たす s, u, v を求める．
template <class T>
std::tuple<Polynomial<T>, Polynomial<T>, Polynomial<T>> Polynomial<T>::extended_gcd(Polynomial<T> f, Polynomial<T> g){
    if(f.deg < g.deg){
        auto tup = Polynomial<T>::extended_gcd(g, f);

        return std::forward_as_tuple(std::get<0>(tup), std::get<2>(tup), std::get<1>(tup));
    }

    Polynomial<T> s = f;
    Polynomial<T> r = g;
    Polynomial<T> u(0, 1);
    Polynomial<T> x(0, 0);
    Polynomial<T> v(0, 0);
    Polynomial<T> y(0, 1);

    while(!r.isZero()){
        auto tup = Polynomial<T>::divide(s, r);
        Polynomial<T> q = std::get<0>(tup);
        Polynomial<T> s1 = r;
        Polynomial<T> r1 = std::get<1>(tup);
        Polynomial<T> u1 = x;
        Polynomial<T> x1 = u - (q * x);
        Polynomial<T> v1 = y;
        Polynomial<T> y1 = v - (q * y);

        s = s1; r = r1; u = u1; x = x1; v = v1; y = y1;
    }

    T lc = s.coeff[s.deg];
    s = s * lc.inv();
    u = u * lc.inv();
    v = v * lc.inv();

    return std::forward_as_tuple(s, u, v);
}

template <class T>
void Polynomial<T>::print() const{
    for(int i = this->deg; i >= 0; --i){
        std::cout << this->coeff[i];
        if(i != 0){
            std::cout << "x^" << i << " + ";
        }
    }
    std::cout << std::endl;
    return;
}

template <class T>
T Polynomial<T>::eval(T x) const{
    T ret(0);
    for(int i = this->deg; i >= 0; --i){
        T c = this->coeff[i];
        for(int j = 1; j < i; ++j){
            c = c * x;
        }
        ret = ret + c;
    }
    return ret;
}

template <class T>
Polynomial<T> Polynomial<T>::derivative() const{
    int deg = this->deg;
    Polynomial<T> r(deg - 1);
    for(int i = 1; i <= this->deg; ++i){
        r.coeff[i - 1] = this->coeff[i] * i;
    }
    return r;
}

template class Polynomial<Number>;
