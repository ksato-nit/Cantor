#pragma once
#include "algorithm"
#include "iostream"

// 抽象的な体の要素を実現する．
class Number {
    static const int CHARA = 31;

public:
    int value;

    Number();
    Number(int value);
    Number operator + (Number y);
    Number operator - (Number y);
    Number operator * (Number y);
    Number operator / (Number y);
    Number inv();

    static Number MINUS_ONE();
};
