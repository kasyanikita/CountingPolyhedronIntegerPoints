#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <gmpxx.h>
#include "../include/ToddPoly_FFT.h"
#include <cmath>

template <typename T>
T random(T range_from, T range_to)
{
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_int_distribution<T> distr(range_from, range_to);
    return distr(generator);
}

template <class TI, class TF>
TF todd3(std::vector<TI> &v, size_t n)
{
    TI quad_pair_sum = 0, triple_sum = 0;
    for (int i = 0; i < n; ++i)
    {
        int rnum = random<int>(0, 50);
        v.push_back(TI(rnum));
    }

    for (int k = 2; k < n; ++k)
    {
        TI pair_sum = 0;
        for (int i = 0; i < k - 1; ++i)
        {
            for (int j = i + 1; j < k; ++j)
            {
                pair_sum = pair_sum + v[i] * v[j];
            }
        }
        triple_sum = triple_sum + v[k] * pair_sum;
    }

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (i != j)
            {
                quad_pair_sum = quad_pair_sum + v[i] * v[i] * v[j];
            }
        }
    }
    return TF(quad_pair_sum) / 24.0 + TF(triple_sum) / 8.0;
}

template <class TI, class TF>
void error(size_t n)
{
    std::vector<mpf_class> error;
    std::vector<mpf_class> errorFFT;
    for (size_t i = 1; i < n; ++i)
    {
        std::vector<TI> v;
        TF ans = todd3<TI, TF>(v, i);
        ToddPoly<TI, TF> todd(3, v);
        ToddPoly_FFT<TI, TF> toddFFT(3, v);
        todd.init();
        toddFFT.init();
        auto t = todd.get_todd();
        auto tFFT = toddFFT.get_todd();
        error.push_back(mpf_class(abs(t[3] - ans), 1000));
        errorFFT.push_back(mpf_class(abs(tFFT[3] - ans), 1000));
    }
    std::ofstream out("data/error.txt");
    for (auto x : error)
    {
        out << mpf_class(x) << std::endl;
    }
    out.close();

    std::ofstream outFFT("data/error_fft.txt");
    for (auto x : errorFFT)
    {
        outFFT << mpf_class(x) << std::endl;
    }
    outFFT.close();
}

int main()
{
    error<mpz_class, mpq_class>(100);
}