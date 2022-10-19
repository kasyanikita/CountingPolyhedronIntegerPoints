#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "../include/GroupElement.h"


TEST_CASE("Create instance") {
    std::vector<GroupNS::intT> comp = {10, 14, 87};
    std::vector<GroupNS::intT> m = {3, 5, 10};
    GroupNS::GroupElement group_elem(m);
    group_elem.assign(comp);

    std::vector<GroupNS::intT> ans = {1, 4, 7};
    REQUIRE(ans == group_elem.get_components());
    REQUIRE(m == group_elem.get_mod());
}


TEST_CASE("Add") {
    std::vector<GroupNS::intT> comp_a = {10, 14, 87};
    std::vector<GroupNS::intT> m = {3, 5, 10};
    GroupNS::GroupElement group_elem_a(m);
    group_elem_a.assign(comp_a);

    std::vector<GroupNS::intT> comp_b = {34, 65, 43};
    GroupNS::GroupElement group_elem_b(m);
    group_elem_b.assign(comp_b);

    auto res = group_elem_a + group_elem_b;

    std::vector<GroupNS::intT> ans = {2, 4, 0};
    REQUIRE(ans == res.get_components());
}


TEST_CASE("Add Equal") {
    std::vector<GroupNS::intT> comp_a = {81, 27, 91};
    std::vector<GroupNS::intT> m = {7, 13, 15};
    GroupNS::GroupElement group_elem_a(m);
    group_elem_a.assign(comp_a);

    std::vector<GroupNS::intT> comp_b = {43, 19, 4};
    GroupNS::GroupElement group_elem_b(m);
    group_elem_b.assign(comp_b);

    group_elem_a += group_elem_b;

    std::vector<GroupNS::intT> ans = {5, 7, 5};
    REQUIRE(ans == group_elem_a.get_components());
}


TEST_CASE("Mod") {
    REQUIRE(GroupNS::modulo(-5, 7) == 2);
    REQUIRE(GroupNS::modulo(-5, -7) == 2);
    REQUIRE(GroupNS::modulo(5, 7) == 5);
    REQUIRE(GroupNS::modulo(5, -7) == 5);
    REQUIRE(GroupNS::modulo(-37, 4) == 3);
    REQUIRE(GroupNS::modulo(-37, -4) == 3);
    REQUIRE(GroupNS::modulo(37, 4) == 1);
    REQUIRE(GroupNS::modulo(37, -4) == 1);
}

TEST_CASE("Subtraction") {
    std::vector<GroupNS::intT> comp_a = {10, 14, 87};
    std::vector<GroupNS::intT> m = {3, 5, 10};
    GroupNS::GroupElement group_elem_a(m);
    group_elem_a.assign(comp_a);

    std::vector<GroupNS::intT> comp_b = {34, 65, 43};
    GroupNS::GroupElement group_elem_b(m);
    group_elem_b.assign(comp_b);

    auto res = group_elem_a - group_elem_b;

    std::vector<GroupNS::intT> ans = {0, 4, 4};
    REQUIRE(ans == res.get_components());
}


TEST_CASE("Scalar multiplication") {
    std::vector<GroupNS::intT> comp = {10, 14, 87};
    std::vector<GroupNS::intT> m = {3, 5, 10};
    GroupNS::GroupElement group_elem(m);
    group_elem.assign(comp);

    GroupNS::intT x = 3;
    std::vector<GroupNS::intT> ans = {0, 2, 1};
    auto res_left_scalar = x * group_elem;
    auto res_right_scalar = group_elem * x; 
    REQUIRE(res_left_scalar.get_components() == ans);
    REQUIRE(res_right_scalar.get_components() == ans);
}


TEST_CASE("Compare false") {
    std::vector<GroupNS::intT> comp_a = {10, 14, 87};
    std::vector<GroupNS::intT> m = {3, 5, 10};
    GroupNS::GroupElement group_elem_a(m);
    group_elem_a.assign(comp_a);

    std::vector<GroupNS::intT> comp_b = {34, 65, 43};
    GroupNS::GroupElement group_elem_b(m);
    group_elem_b.assign(comp_b);

    REQUIRE((group_elem_a < group_elem_b) == false);
}


TEST_CASE("Compare true") {
    std::vector<GroupNS::intT> comp_a = {9, 64, 83};
    std::vector<GroupNS::intT> m = {3, 5, 10};
    GroupNS::GroupElement group_elem_a(m);
    group_elem_a.assign(comp_a);

    std::vector<GroupNS::intT> comp_b = {34, 13, 47};
    GroupNS::GroupElement group_elem_b(m);
    group_elem_b.assign(comp_b);

    REQUIRE((group_elem_a < group_elem_b) == true);
}


TEST_CASE("Get id") {
    std::vector<GroupNS::intT> comp = {10, 14, 87};
    std::vector<GroupNS::intT> m = {3, 5, 10};
    GroupNS::GroupElement group_elem(m);
    group_elem.assign(comp);
    GroupNS::intT ans = 1 + 4 * 3 + 7 * 15;

    REQUIRE(group_elem.get_id() == ans);
}
