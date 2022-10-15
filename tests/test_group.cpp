#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "../include/GroupElement.h"


TEST_CASE("Create instance") {
    std::vector<int> comp = {10, 14, 87};
    std::vector<int> m = {3, 5, 10};
    GroupElement group_elem(comp, m);

    std::vector<int> ans = {1, 4, 7};
    REQUIRE(ans == group_elem.get_components());
    REQUIRE(m == group_elem.get_mod());
}


TEST_CASE("Add") {
    std::vector<int> comp_a = {10, 14, 87};
    std::vector<int> m = {3, 5, 10};
    GroupElement group_elem_a(comp_a, m);

    std::vector<int> comp_b = {34, 65, 43};
    GroupElement group_elem_b(comp_b, m);

    auto res = group_elem_a + group_elem_b;

    std::vector<int> ans = {2, 4, 0};
    REQUIRE(ans == res.get_components());
}


TEST_CASE("Add Equal") {
    std::vector<int> comp_a = {81, 27, 91};
    std::vector<int> m = {7, 13, 15};
    GroupElement group_elem_a(comp_a, m);

    std::vector<int> comp_b = {43, 19, 4};
    GroupElement group_elem_b(comp_b, m);

    group_elem_a += group_elem_b;

    std::vector<int> ans = {5, 7, 5};
    REQUIRE(ans == group_elem_a.get_components());
}
