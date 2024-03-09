cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(CountIntPoints VERSION 1.0)

add_executable(Main src/main.cpp src/counter.cpp
    src/dynamic.cpp src/group_element.cpp
    src/exp_poly.cpp src/tools.cpp src/todd_poly.cpp)

add_executable(Sandbox src/sandbox.cpp src/snf_class.cpp src/counter.cpp
    src/dynamic.cpp src/group_element.cpp src/exp_poly.cpp src/c_vector.cpp
    src/tools.cpp src/todd_poly.cpp src/auxiliary_dynamic.cpp)

include_directories(include)

target_link_libraries(Sandbox PUBLIC gmp gmpxx flint mpfr)
target_link_libraries(Main PUBLIC gmp gmpxx flint mpfr)
