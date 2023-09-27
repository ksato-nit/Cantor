#pragma once

// 抽象的な体の要素を実現する．
class Number {
public:
    int value;
    
    Number(int);
    Number operator + (Number y);
    Number operator - (Number y);
    Number operator * (Number y);
    Number operator / (Number y);
};
