#include "number.hpp"

Number::Number(){
    this->value = 0;
}

Number::Number(int x){
    this->value = x;
}

Number Number::operator + (Number y){
    Number z(y.value + this->value);
    return z;
}

Number Number::operator - (Number y){
    Number z(this->value - y.value);
    return z;
}

Number Number::operator * (Number y){
    Number z(y.value * this->value);
    return z;
}

Number Number::operator / (Number y){
    Number z(this->value / y.value);
    return z;
}
