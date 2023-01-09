#include "tests_defs.hpp"
#include "../include/ExpPoly.h"

TEST_CASE("Multiplication by a number left side")
{
    GroupIP::ExpPoly poly({1, 2, 3, 4}, {1, 4, 10, 7});
    auto res = 2 * poly;
    GroupIP::ExpPoly ans({{1, 2, 3, 4}}, {2, 8, 20, 14});
    REQUIRE(ans.get_poly() == res.get_poly());
}

TEST_CASE("Multiplication by a number right side")
{
    GroupIP::ExpPoly poly({1, 2, 2, 1}, {53, 8, 47, 15});
    auto res = poly * 3;
    GroupIP::ExpPoly ans({1, 2}, {204, 165});
    REQUIRE(ans.get_poly() == res.get_poly());
}

TEST_CASE("Monomial multiplication")
{
    GroupIP::ExpPoly poly({1, 2, 5, 1}, {3, 8, 7, 5});
    auto res = poly.monomial_multiply(7, 4);
    GroupIP::ExpPoly ans({8, 9, 12}, {32, 32, 28});
    REQUIRE(ans.get_poly() == res.get_poly());
}

TEST_CASE("ExpPoly adding")
{
    GroupIP::ExpPoly poly_1({1, 2, 5, 1}, {3, 8, 7, 5});
    GroupIP::ExpPoly poly_2({1, 2, 3, 4}, {1, 4, 10, 7});
    auto res = poly_1 + poly_2;
    GroupIP::ExpPoly ans({1, 2, 3, 4, 5}, {9, 12, 10, 7, 7});
    REQUIRE(ans.get_poly() == res.get_poly());
}