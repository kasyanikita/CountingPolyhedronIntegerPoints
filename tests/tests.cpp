#define CATCH_CONFIG_MAIN
#include <random>
#include <gmpxx.h>
#include "catch.hpp"
#include "../include/ToddPolynomial.h"


template<typename T>
T random(T range_from, T range_to) {
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_int_distribution<T> distr(range_from, range_to);
    return distr(generator);
}


TEST_CASE("Bernulli numbers") {
    Todd<int64_t, mpf_class> t(30, {1, 2});
    t.init();
    auto bernulli = t.get_bernulli();

    REQUIRE(bernulli[0] - 1 < 1e-6);
    REQUIRE(bernulli[1] + 1.0/2 < 1e-6);
    REQUIRE(bernulli[2] - 1.0/6 < 1e-6);
    REQUIRE(bernulli[3] == 0);
    REQUIRE(bernulli[4] + 1.0/30 < 1e-6);
    REQUIRE(bernulli[5] == 0);
    REQUIRE(bernulli[6] - 1.0/42 < 1e-6);
    REQUIRE(bernulli[7] == 0);
    REQUIRE(bernulli[8] + 1.0/30 < 1e-6);
    REQUIRE(bernulli[9] == 0);
    REQUIRE(bernulli[10] - 5.0/66 < 1e-6);
    REQUIRE(bernulli[12] + 691.0/2730 < 1e-6);
    REQUIRE(bernulli[14] - 7.0/6 < 1e-6);
    REQUIRE(bernulli[16] + 3617.0/510 < 1e-6);
    REQUIRE(bernulli[18] - 43867.0/798 < 1e-6);
    REQUIRE(bernulli[20] + 174611.0/330 < 1e-6);
    REQUIRE(bernulli[22] - 854513.0/138 < 1e-6);
    REQUIRE(bernulli[24] + 236364091.0/2730 < 1e-6);
    REQUIRE(bernulli[26] - 8553103.0/6 < 1e-6);
    REQUIRE(bernulli[28] + 23749461029.0/870 < 1e-6);
    REQUIRE(bernulli[30] - 8615841276005.0/14322 < 1e-6);
}


TEST_CASE("Todd Polynomial 0 degree") {
    Todd<int64_t, mpf_class> t1(0, {10});
    t1.init();
    auto todd = t1.get_todd();
    REQUIRE(todd[0] == 1);

    Todd<int64_t, mpf_class> t2(0, {10, 1, 4, 3, 6});
    t2.init();
    todd = t2.get_todd();
    REQUIRE(todd[0] == 1);
}


TEST_CASE("Todd Polynomial 1 degree") {
    Todd<int64_t, mpf_class> t1(1, {65});
    t1.init();
    auto todd = t1.get_todd();
    REQUIRE(todd[0] == 1);
    REQUIRE(todd[1] == 32.5);


    Todd<int64_t, mpf_class> t2(1, {10, 27, 4, 3, 6});
    t2.init();
    todd = t2.get_todd();
    REQUIRE(todd[0] == 1);
    REQUIRE(todd[1] == 25);


    std::vector<int64_t> v;
    int n = 100;
    int64_t sum = 0;
    for (int i = 0; i < n; ++i) {
        int rnum = random<int>(0, 1000);
        sum += rnum;
        v.push_back(rnum);
    }
    Todd<int64_t, mpf_class> t3(1, v);
    t3.init();
    todd = t3.get_todd();
    REQUIRE(todd[0] == 1);
    REQUIRE(todd[1] == static_cast<double>(sum) / 2);
}


TEST_CASE("Todd Polynomial 2 degree") {
    Todd<int64_t, mpf_class> t1(2, {7, 10});
    t1.init();
    auto todd = t1.get_todd();
    REQUIRE(todd[0] == 1);
    REQUIRE(todd[1] == 8.5);
    REQUIRE(todd[2] - 359.0/12 < 1e-6);


    std::vector<int64_t> v;
    int n = 100;
    int64_t sum = 0, quad_sum = 0, mul_sum = 0;
    for (int i = 0; i < n; ++i) {
        int rnum = random<int>(0, 1000);
        sum += rnum;
        quad_sum += rnum * rnum;
        v.push_back(rnum);
    }
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            mul_sum += v[i] * v[j];
        }
    }
    Todd<int64_t, mpf_class> t3(1, v);
    t3.init();
    todd = t3.get_todd();
    REQUIRE(todd[0] == 1);
    REQUIRE(todd[1] == sum / 2.0);
    REQUIRE(todd[2] - (quad_sum / 12.0 + mul_sum / 4.0) < 1e-6);
}

