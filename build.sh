#!/bin/bash

# build ディレクトリがすでにあれば削除する
if [ -e ./build ]; then
    rm -rf ./build
fi

mkdir ./build

cd ./build

cmake ..

make
