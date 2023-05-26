#pragma once

#include <random>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

const auto tests_seed = std::random_device()();
std::mt19937 tests_rand_gen(tests_seed);

template <typename T>
T tests_rand(T range_from, T range_to)
{
    std::uniform_int_distribution<T> udistr(range_from, range_to);
    return udistr(tests_rand_gen);
}  // delete