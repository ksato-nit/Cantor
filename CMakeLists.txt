# 要求する CMake のバージョン
cmake_minimum_required(VERSION 3.5.0)

# C++ 14 以上を要求
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Mac で GCC の使用を強制
if(CMAKE_SYSTEM_NAME MATCHES "Sonoma")
    set(CMAKE_C_COMPILER "/opt/homebrew/bin/gcc" CACHE string "clang compiler" FORCE)
    set(CMAKE_CXX_COMPILER "/opt/homebrew/bin/g++" CACHE string "clang++ compiler" FORCE)
endif()

# プロジェクト名
project(Cantor)

# ディレクトリ設定
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(INC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/inc)
set(TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests)

# GoogleTest の設定
include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)

find_package(Boost REQUIRED)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# コンパイルオプション設定
set(CMAKE_CXX_FLAGS "-O3 -std=c++14 -Wall -lgmpxx -lgmp -pg")

# ワイルドカード指定の設定
file(GLOB SRC ${SRC_DIR}/*.cpp)
file(GLOB TEST ${TEST_DIR}/*.cpp)

include_directories(
    ${INC_DIR}
    ${Boost_INCLUDE_DIRS}
)

add_executable(
    Cantor
    ${SRC}
)

target_link_libraries(
    Cantor
    gmp
    gmpxx
)

# テストの設定
enable_testing()

add_executable(
    CantorTest
    ${SRC_DIR}/number.cpp ${SRC_DIR}/polynomial.cpp ${SRC_DIR}/mumford.cpp ${SRC_DIR}/mumford_projective.cpp ${SRC_DIR}/mumford_weighted_projective.cpp ${SRC_DIR}/extended_number.cpp
    ${TEST}
)
target_link_libraries(
    CantorTest
    GTest::gtest_main
    gmp
    gmpxx
)
target_include_directories(
    CantorTest
    PUBLIC
    ${INC_DIR}
)

include(GoogleTest)
gtest_discover_tests(CantorTest)
