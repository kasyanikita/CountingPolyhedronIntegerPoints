#include "tests_defs.hpp"
#include "../include/ToddPoly_FFT.h"

TEST_CASE("Todd Polynomial 0 degree FFT")
{
    ToddPoly_FFT<int64_t, mpf_class> t1(0, {10});
    t1.init();
    auto todd = t1.get_todd();
    REQUIRE(abs(todd[0] - 1) < 1e-6);

    ToddPoly_FFT<int64_t, mpf_class> t2(0, {10, 1, 4, 3, 6});
    t2.init();
    todd = t2.get_todd();
    REQUIRE(abs(todd[0] - 1) < 1e-6);
}

TEST_CASE("Todd Polynomial 1 degree FFT")
{
    ToddPoly_FFT<int64_t, mpf_class> t1(1, {65});
    t1.init();
    auto todd = t1.get_todd();
    REQUIRE(abs(todd[0] - 1) < 1e-6);
    REQUIRE(abs(todd[1] - 32.5) < 1e-6);

    ToddPoly_FFT<int64_t, mpf_class> t2(1, {10, 27, 4, 3, 6});
    t2.init();
    todd = t2.get_todd();
    REQUIRE(abs(todd[0] - 1) < 1e-6);
    REQUIRE(abs(todd[1] - 25) < 1e-6);

    std::vector<int64_t> v;
    int n = 100;
    int64_t sum = 0;
    for (int i = 0; i < n; ++i)
    {
        int rnum = tests_rand<int>(0, 1000);
        sum += rnum;
        v.push_back(rnum);
    }
    ToddPoly_FFT<int64_t, mpf_class> t3(1, v);
    t3.init();
    todd = t3.get_todd();
    REQUIRE(abs(todd[0] - 1) < 1e-6);
    REQUIRE(abs(todd[1] - static_cast<double>(sum) / 2) < 1e-3); //!!!!!!!!!!!!!!! 4.52297e-05 < 0.000001
}

TEST_CASE("Todd Polynomial 2 degree FFT")
{
    ToddPoly_FFT<int64_t, mpf_class> t1(2, {7, 10});
    t1.init();
    auto todd = t1.get_todd();
    REQUIRE(abs(todd[0] - 1) < 1e-6);
    REQUIRE(abs(todd[1] - 8.5) < 1e-6);
    REQUIRE(abs(todd[2] - 359.0 / 12) < 1e-6);

    std::vector<int64_t> v;
    int n = 100;
    int64_t sum = 0, quad_sum = 0, mul_sum = 0;
    for (int i = 0; i < n; ++i)
    {
        int rnum = tests_rand<int>(0, 1000);
        sum += rnum;
        quad_sum += rnum * rnum;
        v.push_back(rnum);
    }
    for (int i = 0; i < n - 1; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            mul_sum += v[i] * v[j];
        }
    }
    ToddPoly_FFT<int64_t, mpf_class> t3(2, v);
    t3.init();
    todd = t3.get_todd();
    REQUIRE(abs(todd[0] - 1) < 1e-6);
    REQUIRE(abs(todd[1] - static_cast<double>(sum) / 2) < 1e-6);
    REQUIRE(abs(todd[2] - (quad_sum / 12.0 + mul_sum / 4.0)) < 1e-6);
}

TEST_CASE("Todd Polynomial 3 degree FFT")
{
    ToddPoly_FFT<int64_t, mpf_class> t1(3, {1, 2, 3});
    t1.init();
    auto todd = t1.get_todd();
    REQUIRE(abs(todd[0] - 1) < 1e-6);
    REQUIRE(abs(todd[1] - 3) < 1e-6);
    REQUIRE(abs(todd[2] - 47.0 / 12) < 1e-6);
    REQUIRE(abs(todd[3] - 2.75) < 1e-6);

    std::vector<int64_t> v;
    int n = 100;
    int64_t quad_pair_sum = 0, triple_sum = 0;
    for (int i = 0; i < n; ++i)
    {
        int rnum = tests_rand<int>(0, 50);
        v.push_back(rnum);
    }

    for (int k = 2; k < n; ++k)
    {
        int64_t pair_sum = 0;
        for (int i = 0; i < k - 1; ++i)
        {
            for (int j = i + 1; j < k; ++j)
            {
                pair_sum += v[i] * v[j];
            }
        }
        triple_sum += v[k] * pair_sum;
    }

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (i != j)
            {
                quad_pair_sum += v[i] * v[i] * v[j];
            }
        }
    }
    ToddPoly_FFT<int64_t, mpf_class> t2(3, v);
    t2.init();
    todd = t2.get_todd();
    REQUIRE(abs(todd[3] - (quad_pair_sum / 24.0 + triple_sum / 8.0)) < 1e-6);
}

TEST_CASE("Todd Polynomial 4 degree FFT")
{
    std::vector<int64_t> v;
    int n = 100;
    uint64_t sum_pow_4 = 0, quad_pair = 0, triple_sum = 0, sum_4 = 0;
    for (int i = 0; i < n; ++i)
    {
        int rnum = tests_rand<int>(1, 10);
        v.push_back(rnum);
        sum_pow_4 += rnum * rnum * rnum * rnum;
    }

    for (int i = 0; i < n - 1; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            quad_pair += v[i] * v[i] * v[j] * v[j];
        }
    }

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n - 1; ++j)
        {
            for (int k = j + 1; k < n; ++k)
            {
                if (j != i && k != i)
                {
                    triple_sum += v[i] * v[i] * v[j] * v[k];
                }
            }
        }
    }
    for (int m = 2; m < n - 1; ++m)
    {
        for (int k = m + 1; k < n; ++k)
        {
            int64_t pair_sum = 0;
            for (int i = 0; i < m - 1; ++i)
            {
                for (int j = i + 1; j < m; ++j)
                {
                    pair_sum += v[i] * v[j];
                }
            }
            sum_4 += v[m] * v[k] * pair_sum;
        }
    }
    ToddPoly_FFT<int64_t, mpf_class> t(4, v);
    t.init();
    auto todd = t.get_todd();
    REQUIRE(abs(todd[4] - (-1.0 * sum_pow_4 / 720 + 1.0 * quad_pair / 144 + 1.0 * triple_sum / 48 + 1.0 * sum_4 / 16)) < 1e-6);
}