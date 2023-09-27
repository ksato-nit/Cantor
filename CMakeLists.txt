# 要求する CMake のバージョン
cmake_minimum_required(VERSION 3.5.0)

# C++ 11 以上を要求
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# プロジェクト名
project(Cantor)

# ディレクトリ設定
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(INC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/inc)

# コンパイルオプション設定
set(CMAKE_CXX_FLAGS "-O2 -std=c++11 -Wall")

# ワイルドカード指定の設定
file(GLOB SRC ${SRC_DIR}/*.cpp)

include_directories(
    ${MPI_INCLUDE_PATH}
    ${INC_DIR}
)

add_executable(
    Cantor
    ${SRC}
)
