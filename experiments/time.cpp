#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include <gmpxx.h>
#include "../include/ToddPoly_FFT.h"

template <typename T>
T random(T range_from, T range_to)
{
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_int_distribution<T> distr(range_from, range_to);
    return distr(generator);
}

template <class TI, class TF>
void time_fixed_input(size_t m, size_t n)
{
    std::vector<TI> vfft;
    for (int i = 0; i < n; ++i)
    {
        int rnum = random<int>(0, 100);
        vfft.push_back(rnum);
    }
    std::vector<TF> times;
    for (int i = 1; i <= m; ++i)
    {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        ToddPoly_FFT<TI, TF> tfft(i, vfft);
        tfft.init();
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        times.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9);
    }
    std::ofstream fftout("data/time_fft_fixed_input.txt");
    for (auto x : times)
    {
        fftout << mpf_class(x) << std::endl;
    }
    fftout.close();

    times.clear();
    std::vector<TI> v;
    for (int i = 0; i < n; ++i)
    {
        int rnum = random<int>(0, 100);
        v.push_back(rnum);
    }
    for (int i = 1; i <= m; ++i)
    {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        Todd<TI, TF> t(i, v);
        t.init();
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        times.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9);
    }
    std::ofstream out("data/time_fixed_input.txt");
    for (auto x : times)
    {
        out << mpf_class(x) << std::endl;
    }
    out.close();
}

template <class TI, class TF>
void time_fixed_degree(size_t m, size_t n)
{
    std::vector<TF> times;
    std::vector<TI> v = {1};
    for (int i = 1; i <= n; ++i)
    {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        Todd<TI, TF> t(m, v);
        t.init();
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        times.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9);
        v.push_back(random<int>(1, 100));
    }
    std::ofstream out("data/time_fixed_degree.txt");
    for (auto x : times)
    {
        out << mpf_class(x) << std::endl;
    }
    out.close();
}

int main()
{
    time_fixed_input<mpz_class, mpq_class>(150, 10);
    time_fixed_degree<mpz_class, mpq_class>(10, 100);
}