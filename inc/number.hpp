#pragma once

// 抽象的な体の要素を実現する．
class Number {
public:
    double value;

    Number();
    Number(double);
    Number operator + (Number y);
    Number operator - (Number y);
    Number operator * (Number y);
    Number operator / (Number y);
};
