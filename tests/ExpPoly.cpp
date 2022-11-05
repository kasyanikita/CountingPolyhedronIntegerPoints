#include "tests_defs.hpp"
#include "../include/ExpPoly.h"

TEST_CASE("Multiplication by a number left side")
{
    GroupIP::ExpPoly poly({1, 2, 3, 4});
    auto res = 2 * poly;
    std::vector<double> ans = {2, 2, 2, 2};
    REQUIRE(ans == res.get_base());
}

TEST_CASE("Multiplication by a number right side")
{
    GroupIP::ExpPoly poly({1, 2, 3, 4});
    auto res = poly * 12;
    std::vector<double> ans = {12, 12, 12, 12};
    REQUIRE(ans == res.get_base());
}