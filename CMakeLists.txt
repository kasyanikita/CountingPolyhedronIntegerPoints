cmake_minimum_required(VERSION 3.16)

project(CountIntPoints VERSION 1.0 LANGUAGES CXX)

add_compile_options(-I ../include)
add_executable(CountIntPoints src/main.cpp src/ToddPolynomial.cpp)