#include "number.hpp"

Number::Number(){
    this->value = 0;
}

Number::Number(int x){
    this->value = x % CHARA;
}

Number Number::operator + (Number y){
    int value = (y.value + this->value) % CHARA;
    Number z(value);
    return z;
}

Number Number::operator - (Number y){
    int value = (y.value - this->value) % CHARA;
    Number z(value);
    return z;
}

Number Number::operator * (Number y){
    int value = (y.value * this->value) % CHARA;
    Number z(value);
    return z;
}

Number Number::operator / (Number y){
    Number yinv = y.inv();
    int value = (this->value * yinv.value) % CHARA;
    Number z(value);
    return z;
}

Number Number::inv(){
    if(this->value < 0){
        Number mx(this->value * -1);
        return mx.inv();
    }

    int s = this->value, t = CHARA;
    int x = 1, u = 0;
    int k;
    while(t > 0) {
        k = s / t;
        s = s - (k * t);
        std::swap(s, t);
        x = x - (k * u);
        std::swap(x, u);
    }
    x = x % CHARA;
    if(x < 0){
        x = x + CHARA;
    }
    return x;
}
