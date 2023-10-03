#include "polynomial.hpp"

Polynomial::Polynomial(){
    this->deg = 0;
    this->coeff.resize(1);
    this->coeff[0] = 0;
}

Polynomial::Polynomial(int deg){
    this->deg = deg;
    this->coeff.resize(deg + 1);
}

Polynomial::Polynomial(int deg, int* coeff){
    this->deg = deg;
    this->coeff.resize(deg + 1);
    for(int i = 0; i <= deg; ++i){
        Number c(coeff[i]);
        this->coeff[i] = c;
    }
}

Polynomial::Polynomial(int deg, int c0){
    this->deg = deg;
    this->coeff.resize(deg + 1);
    Number c(c0);
    this->coeff[0] = c;
}

Polynomial::Polynomial(int deg, Number c0){
    this->deg = deg;
    this->coeff.resize(deg + 1);
    this->coeff[0] = c0;
}

Polynomial::Polynomial(int deg, Number c0, Number c1){
    this->deg = deg;
    this->coeff.resize(deg + 1);
    this->coeff[0] = c0;
    this->coeff[1] = c1;
}

void Polynomial::resize(int deg){
    this->deg = deg;
    this->coeff.resize(deg + 1);
}

void Polynomial::resize(std::vector<Number> coeff){
    this->deg = deg;
    this->coeff = coeff;
}

// 次数が整合的であるようにする (高次の係数で 0 のものがあれば切り落としていく)．
void Polynomial::normalize(){
    int deg = this->deg;
    for(auto i = this->coeff.rbegin(), e = this->coeff.rend(); i != e; ++i){
        if(i->value == 0){
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

bool Polynomial::isZero(){
    for(auto n : this->coeff){
        if(n.value != 0){
            return false;
        }
    }
    return true;
}

Polynomial Polynomial::operator + (Polynomial q){
    // TODO : もう少しスマートな実装を考える．
    int deg = std::max(this->deg, q.deg);
    Polynomial r(deg);
    for(int i = 0; i <= deg; ++i){
        Number a = (this->deg >= i) ? this->coeff[i] : 0;
        Number b = (q.deg >= i) ? q.coeff[i] : 0;
        r.coeff[i] = a + b;
    }
    //r.normalize();
    return r;
}

Polynomial Polynomial::operator - (Polynomial q){
    int deg = std::max(this->deg, q.deg);
    Polynomial r(deg);
    for(int i = 0; i <= deg; ++i){
        Number a = (this->deg >= i) ? this->coeff[i] : 0;
        Number b = (q.deg >= i) ? q.coeff[i] : 0;
        r.coeff[i] = a - b;
    }
    //r.normalize();
    return r;
}

Polynomial Polynomial::operator * (Polynomial q){
    int deg = this->deg + q.deg;
    Polynomial r(deg);
    for(int i = 0; i <= this->deg; ++i){
        for(int j = 0; j <= q.deg; ++j){
            Number a = (this->deg >= i) ? this->coeff[i] : 0;
            Number b = (q.deg >= j) ? q.coeff[j] : 0;
            r.coeff[i + j] = r.coeff[i + j] + a * b;
        }
    }
    //r.normalize();
    return r;
}

Polynomial Polynomial::operator * (Number k){
    int deg = this->deg;
    Polynomial r(deg);
    for(int i = 0; i <= this->deg; ++i){
        r.coeff[i] = this->coeff[i] * k;
    }
    return r;
}

Polynomial Polynomial::operator / (Polynomial g){
    auto tup = Polynomial::divide(*this, g);

    Polynomial q = std::get<0>(tup);
    Polynomial r = std::get<1>(tup);

    return q;
}

Polynomial Polynomial::operator % (Polynomial g){
    if(this->deg < g.deg){
        return *this;
    }

    auto tup = Polynomial::divide(*this, g);

    Polynomial q = std::get<0>(tup);
    Polynomial r = std::get<1>(tup);

    return r;
}

// f = qg + r を満たす q, r を求める．
// deg g <= deg f を仮定．
std::tuple<Polynomial, Polynomial> Polynomial::divide(Polynomial f, Polynomial g){
    int deg_f = f.deg;
    int deg_g = g.deg;
    Polynomial q(deg_f - deg_g);
    Polynomial r(deg_f);

    // TODO : インデックスをシフトしてもう少しきれいに書けそう．
    for(int j = deg_f - deg_g; j >= 0; --j){
        q.coeff[j] = f.coeff[j + deg_g] / g.coeff[deg_g];
        Polynomial mono(j);
        mono.coeff[j] = q.coeff[j];
        r = f - (g * mono);
        f = r;
    }

    r.normalize();

    return std::forward_as_tuple(q, r);
}


// s = gcd(f, g), fu + gv = s を満たす s, u, v を求める．
std::tuple<Polynomial, Polynomial, Polynomial> Polynomial::extended_gcd(Polynomial f, Polynomial g){
    if(f.deg < g.deg){
        auto tup = Polynomial::extended_gcd(g, f);

        return std::forward_as_tuple(std::get<0>(tup), std::get<2>(tup), std::get<1>(tup));
    }

    Polynomial s = f;
    Polynomial r = g;
    Polynomial u(0, 1);
    Polynomial x(0, 0);
    Polynomial v(0, 0);
    Polynomial y(0, 1);

    while(!r.isZero()){
        auto tup = Polynomial::divide(s, r);
        Polynomial q = std::get<0>(tup);
        Polynomial s1 = r;
        Polynomial r1 = std::get<1>(tup);
        Polynomial u1 = x;
        Polynomial x1 = u - (q * x);
        Polynomial v1 = y;
        Polynomial y1 = v - (q * y);

        s = s1; r = r1; u = u1; x = x1; v = v1; y = y1;
    }

    Number lc = s.coeff[s.deg];
    s = s * lc.inv();
    u = u * lc.inv();
    v = v * lc.inv();

    return std::forward_as_tuple(s, u, v);
}

void Polynomial::print(){
    for(int i = this->deg; i >= 0; --i){
        std::cout << this->coeff[i].value;
        if(i != 0){
            std::cout << "x^" << i << " + ";
        }
    }
    std::cout << std::endl;
}
