cmake_minimum_required(VERSION 3.16)

project(CountIntPoints)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
include_directories(CountIntPoints "./include")

add_executable(CountIntPoints src/main.cpp)
target_link_libraries(CountIntPoints gmp gmpxx flint mpfr)