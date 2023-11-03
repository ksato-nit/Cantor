# 要求する CMake のバージョン
cmake_minimum_required(VERSION 3.5.0)

# C++ 14 以上を要求
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# プロジェクト名
project(Cantor)

# ディレクトリ設定
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(INC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/inc)

include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# コンパイルオプション設定
set(CMAKE_CXX_FLAGS "-O2 -std=c++14 -Wall")

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
