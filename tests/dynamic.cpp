#include "tests_defs.hpp"
#include "../include/Dynamic.h"
using namespace GroupIP;

TEST_CASE("Get matrix column")
{
    std::vector<std::vector<int_t>> M = {
        {1, 2, 3, 4},
        {5, 6, 7, 8},
        {9, 10, 11, 12},
        {13, 14, 15, 16}
    };
    std::vector<int_t> ans = {3, 7, 11, 15};
    auto res = get_mat_col(M, 2);
    REQUIRE(ans == res);
}

TEST_CASE("Matrix vector multiplication")
{
    std::vector<std::vector<int_t>> M = {
        {1, 2, 3, 4},
        {5, 6, 7, 8},
        {9, 10, 11, 12},
        {13, 14, 15, 16}
    };
    std::vector<int_t> v = {1, -1, 5, -3};
    auto res = mat_vec_mult(M, v);
    std::vector<int_t> ans = {2, 10, 18, 26};
    REQUIRE(res == ans);
}

TEST_CASE("Calculation g vector")
{
    std::vector<std::vector<int_t>> P = {
        {1, 0, 2, 4},
        {7, 1, 0, 8},
        {2, 7, 3, 9},
        {4, 1, 8, 0}
    };
    std::vector<int_t> b = {6, 5, 10, 3};
    std::vector<int_t> S = {1, 1, 2, 8};
    auto res = calc_g(P, b, S);

    // check size of g vec
    REQUIRE(res.size() == 5);

    // check mod vector
    REQUIRE(res[0].get_mod() == S);
    REQUIRE(res[1].get_mod() == S);
    REQUIRE(res[2].get_mod() == S);
    REQUIRE(res[3].get_mod() == S);
    REQUIRE(res[4].get_mod() == S);

    // check components
    REQUIRE(res[0].get_components() == std::vector<int_t>({0, 0, 0, 5}));
    REQUIRE(res[1].get_components() == std::vector<int_t>({0, 0, 0, 4}));
    REQUIRE(res[2].get_components() == std::vector<int_t>({0, 0, 1, 1}));
    REQUIRE(res[3].get_components() == std::vector<int_t>({0, 0, 1, 0}));
    REQUIRE(res[4].get_components() == std::vector<int_t>({0, 0, 1, 0}));
}

TEST_CASE("Calculation r vector") {
    std::vector<std::vector<int_t>> P = {
        {1, 0, 2, 4},
        {7, 1, 0, 8},
        {2, 7, 3, 9},
        {4, 1, 8, 0}};
    std::vector<int_t> b = {6, 5, 10, 3};
    std::vector<int_t> S = {1, 1, 2, 8};
    auto g = calc_g(P, b, S);
    auto res = calc_r(g);
    REQUIRE(res == std::vector<uint_t>({8, 2, 8, 2, 2}));
}

TEST_CASE("Calculation s vector")
{
    Dynamic d;
    std::vector<GroupElement> g(2, GroupElement({1, 1, 2, 8}));
    g[0].assign({0, 0, 0, 4});
    g[1].assign({0, 0, 0, 5});
    std::vector<uint_t> r({2, 8});
    auto s = calc_s(g, r);
    REQUIRE(s == std::vector<uint_t>({4, 1}));
}

TEST_CASE("Calculation h vector")
{
    std::vector<std::vector<int_t>> A = {
        {1, 9, -5},
        {0, 6, 0},
        {0, 3, 6}};
    auto res = calc_h(A);
    std::vector<std::vector<int_t>> ans = {
        {36, 0, 0},
        {-69, 6, -3},
        {30, 0, 6}};
    REQUIRE(ans == res);
}