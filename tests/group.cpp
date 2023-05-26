#include "tests_defs.hpp"
#include "../include/GroupElement.h"
#include "Dynamic.h"

using namespace GroupIP;

TEST_CASE("Create instance")
{
    std::vector<GroupIP::int_t> comp = {10, 14, 87};
    std::vector<GroupIP::int_t> m = {3, 5, 10};
    GroupElement group_elem(m);
    group_elem.assign(comp);

    std::vector<GroupIP::int_t> ans = {1, 4, 7};
    REQUIRE(ans == group_elem.get_components());
    REQUIRE(m == group_elem.get_mod());
}

TEST_CASE("Add")
{
    std::vector<GroupIP::int_t> comp_a = {10, 14, 87};
    std::vector<GroupIP::int_t> m = {3, 5, 10};
    GroupElement group_elem_a(m);
    group_elem_a.assign(comp_a);

    std::vector<GroupIP::int_t> comp_b = {34, 65, 43};
    GroupElement group_elem_b(m);
    group_elem_b.assign(comp_b);

    auto res = group_elem_a + group_elem_b;

    std::vector<GroupIP::int_t> ans = {2, 4, 0};
    REQUIRE(ans == res.get_components());
}

TEST_CASE("Add Equal")
{
    std::vector<GroupIP::int_t> comp_a = {81, 27, 91};
    std::vector<GroupIP::int_t> m = {7, 13, 15};
    GroupElement group_elem_a(m);
    group_elem_a.assign(comp_a);

    std::vector<GroupIP::int_t> comp_b = {43, 19, 4};
    GroupElement group_elem_b(m);
    group_elem_b.assign(comp_b);

    group_elem_a += group_elem_b;

    std::vector<GroupIP::int_t> ans = {5, 7, 5};
    REQUIRE(ans == group_elem_a.get_components());
}

TEST_CASE("Mod")
{
    REQUIRE(GroupIP::modulo(-5, 7) == 2);
    // REQUIRE(GroupIP::modulo(-5, -7) == 2);
    REQUIRE(GroupIP::modulo(5, 7) == 5);
    // REQUIRE(GroupIP::modulo(5, -7) == 5);
    REQUIRE(GroupIP::modulo(-37, 4) == 3);
    // REQUIRE(GroupIP::modulo(-37, -4) == 3);
    REQUIRE(GroupIP::modulo(37, 4) == 1);
    // REQUIRE(GroupIP::modulo(37, -4) == 1);
}

TEST_CASE("Subtraction")
{
    std::vector<GroupIP::int_t> comp_a = {10, 14, 87};
    std::vector<GroupIP::int_t> m = {3, 5, 10};
    GroupElement group_elem_a(m);
    group_elem_a.assign(comp_a);

    std::vector<GroupIP::int_t> comp_b = {34, 65, 43};
    GroupElement group_elem_b(m);
    group_elem_b.assign(comp_b);

    auto res = group_elem_a - group_elem_b;

    std::vector<GroupIP::int_t> ans = {0, 4, 4};
    REQUIRE(ans == res.get_components());
}

TEST_CASE("Scalar multiplication")
{
    std::vector<GroupIP::int_t> comp = {10, 14, 87};
    std::vector<GroupIP::int_t> m = {3, 5, 10};
    GroupElement group_elem(m);
    group_elem.assign(comp);

    GroupIP::int_t x = 3;
    std::vector<GroupIP::int_t> ans = {0, 2, 1};
    auto res_left_scalar = x * group_elem;
    auto res_right_scalar = group_elem * x;
    REQUIRE(res_left_scalar.get_components() == ans);
    REQUIRE(res_right_scalar.get_components() == ans);
}

TEST_CASE("Compare false")
{
    std::vector<GroupIP::int_t> comp_a = {10, 14, 87};
    std::vector<GroupIP::int_t> m = {3, 5, 10};
    GroupElement group_elem_a(m);
    group_elem_a.assign(comp_a);

    std::vector<GroupIP::int_t> comp_b = {34, 65, 43};
    GroupElement group_elem_b(m);
    group_elem_b.assign(comp_b);

    REQUIRE((group_elem_a < group_elem_b) == false);
}

TEST_CASE("Compare true")
{
    std::vector<GroupIP::int_t> comp_a = {9, 64, 83};
    std::vector<GroupIP::int_t> m = {3, 5, 10};
    GroupElement group_elem_a(m);
    group_elem_a.assign(comp_a);

    std::vector<GroupIP::int_t> comp_b = {34, 13, 47};
    GroupElement group_elem_b(m);
    group_elem_b.assign(comp_b);

    REQUIRE((group_elem_a < group_elem_b) == true);
}

TEST_CASE("Get id")
{
    std::vector<GroupIP::int_t> comp = {10, 14, 87};
    std::vector<GroupIP::int_t> m = {3, 5, 10};
    GroupElement group_elem(m);
    group_elem.assign(comp);
    GroupIP::int_t ans = 1 + 4 * 3 + 7 * 15;

    REQUIRE(group_elem.get_idx() == ans);
}

TEST_CASE("Get group element by index") {
    using namespace GroupIP;
    std::vector<int_t> S = {1, 2, 816};
    auto res = get_group_element_by_index(537, S);
    std::vector<int_t> ans({0, 1, 268});
    REQUIRE(res.get_components() == ans);
}
