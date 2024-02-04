# Cantor
[![build](https://github.com/ksato-nit/Cantor/actions/workflows/build.yml/badge.svg)](https://github.com/ksato-nit/Cantor/actions/workflows/build.yml)

超楕円曲線上の加算アルゴリズムの C++ による実装．

2023 年度卒業研究の一環として作った．

# Branches
- main : ハイレベル実装，素体
- mpz_t : 低レベル実装，素体
  - 整数演算のオーバーロードを使わない実装
- extended_field : 低レベル実装，2 次拡大体
  - 卒業論文に記載の数値実験に使用

# Dependencies
- GCC 3.11
- CMake 3.5
- Boost C++
- GMP

# Build
- ``` mkdir build && cd build && cmake .. && make```
または，
- ```./build.sh```
