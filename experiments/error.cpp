#include <iostream>
#include <vector>
#include <random>
// #include <chrono>
#include <fstream>
#include <gmpxx.h>
#include "../include/ToddFFT.h"
#include <cmath>

template<typename T>
T random(T range_from, T range_to) {
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_int_distribution<T> distr(range_from, range_to);
    return distr(generator);
}

mpq_class todd3(std::vector<mpz_class>& v, size_t n) {
    mpz_class quad_pair_sum = 0, triple_sum = 0;
    for (int i = 0; i < n; ++i) {
        int rnum = random<int>(0, 50);
        v.push_back(mpz_class(rnum));
    }
    
    for (int k = 2; k < n; ++k) {
        mpz_class pair_sum = 0;
        for (int i = 0; i < k - 1; ++i) {
            for (int j = i + 1; j < k; ++j) {
                pair_sum = pair_sum + v[i] * v[j];
            }
        }
        triple_sum = triple_sum + v[k] * pair_sum;
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                quad_pair_sum = quad_pair_sum + v[i] * v[i] * v[j];
            }
        }
    }
    return mpq_class(quad_pair_sum.get_d() / 24.0 + triple_sum.get_d() / 8.0);
}

void error(size_t n) {
    std::vector<mpf_class> error;
    std::vector<mpf_class> errorFFT;
    for (size_t i = 1; i < n; ++i) {
        std::vector<mpz_class> v;
        mpq_class ans = todd3(v, i);
        Todd<mpz_class, mpq_class> todd(3, v);
        ToddFFT<mpz_class, mpq_class> toddFFT(3, v);
        todd.init();
        toddFFT.init();
        auto t = todd.get_todd();
        auto tFFT = toddFFT.get_todd();
        error.push_back(mpf_class(abs(t[3] - ans), 1000));
        errorFFT.push_back(mpf_class(abs(tFFT[3] - ans), 1000));
    }
    std::ofstream out("../data/error.txt");
    for (auto x : error) {
        out << mpf_class(x) << std::endl;
    }
    out.close();

    std::ofstream outFFT("../data/error_fft.txt");
    for (auto x : errorFFT) {
        outFFT << mpf_class(x) << std::endl;
    }
    outFFT.close();
}

int main() {
    error(200);
}