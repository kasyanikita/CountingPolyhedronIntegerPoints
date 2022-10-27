cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(CountIntPoints VERSION 1.0)

add_executable(CountIntPoints src/main.cpp)
add_executable(CountIntPoints_Tests tests/test_all.cpp)

target_include_directories(CountIntPoints PUBLIC "./include")
target_link_libraries(CountIntPoints PUBLIC gmp gmpxx flint mpfr)

target_include_directories(CountIntPoints_Tests PUBLIC "./include")
target_link_libraries(CountIntPoints_Tests PUBLIC gmp gmpxx flint mpfr)