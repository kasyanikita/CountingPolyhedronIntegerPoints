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

double todd3(std::vector<int64_t>& v, size_t n) {
    int64_t quad_pair_sum = 0, triple_sum = 0;
    for (int i = 0; i < n; ++i) {
        int rnum = random<int>(0, 50);
        v.push_back(rnum);
    }
    
    for (int k = 2; k < n; ++k) {
        int64_t pair_sum = 0;
        for (int i = 0; i < k - 1; ++i) {
            for (int j = i + 1; j < k; ++j) {
                pair_sum += v[i] * v[j];
            }
        }
        triple_sum += v[k] * pair_sum;
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                quad_pair_sum += v[i] * v[i] * v[j];
            }
        }
    }
    return quad_pair_sum / 24.0 + triple_sum / 8.0;
}

void error(size_t n) {
    std::vector<double> error;
    std::vector<double> errorFFT;
    for (size_t i = 1; i < n; ++i) {
        std::vector<int64_t> v;
        double ans = todd3(v, i);
        Todd<int64_t, double> todd(3, v);
        ToddFFT<int64_t, double> toddFFT(3, v);
        todd.init();
        toddFFT.init();
        auto t = todd.get_todd();
        auto tFFT = toddFFT.get_todd();
        error.push_back(std::abs(t[3] - ans));
        errorFFT.push_back(std::abs(tFFT[3] - ans));
    }
    std::ofstream out("../data/error.txt");
    for (auto x : error) {
        out << x << std::endl;
    }
    out.close();

    std::ofstream outFFT("../data/error_fft.txt");
    for (auto x : errorFFT) {
        outFFT << x << std::endl;
    }
    outFFT.close();
}

int main() {
    error(100);
}