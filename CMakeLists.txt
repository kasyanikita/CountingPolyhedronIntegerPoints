cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(CountIntPoints VERSION 1.0)

add_executable(CountIntPoints src/main.cpp)
add_executable(CountIntPoints_Tests tests/test_all.cpp)
add_executable(SNF_stress_test tests/stress_snf.cpp)
add_executable(ConstantTailor src/constant_in_tailor.cpp)
add_executable(Stress_tests tests/stress.cpp)
add_executable(Time experiments/time.cpp)


target_include_directories(ConstantTailor PUBLIC "./include")
target_link_libraries(ConstantTailor PUBLIC gmp gmpxx flint mpfr)

target_include_directories(CountIntPoints PUBLIC "./include")
target_link_libraries(CountIntPoints PUBLIC gmp gmpxx flint mpfr)

target_include_directories(CountIntPoints_Tests PUBLIC "./include")
target_link_libraries(CountIntPoints_Tests PUBLIC gmp gmpxx flint mpfr)

target_include_directories(SNF_stress_test PUBLIC "./include")
target_link_libraries(SNF_stress_test PUBLIC gmp gmpxx flint mpfr)

target_include_directories(Stress_tests PUBLIC "./include")
target_link_libraries(Stress_tests PUBLIC gmp gmpxx flint mpfr)

target_include_directories(Time PUBLIC "./include")
target_link_libraries(Time PUBLIC gmp gmpxx flint mpfr)