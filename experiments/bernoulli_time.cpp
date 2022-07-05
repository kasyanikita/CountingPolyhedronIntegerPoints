#include <fstream>
#include <iostream>
#include <chrono>
#include <vector>

std::vector<int64_t> calc_pascal(size_t n) {
    // Get Pascal's triangle n-th row
    std::vector<int64_t> res(n + 1, 1);
    for (size_t i = 0; i < n; ++i) {
        res[i + 1] = res[i] * (n - i) / (i + 1);
    }
    return res;
}


std::vector<double> calc_bernoulli(size_t n) {
    // Calculate bernoulli numbers from 0 to m
    std::vector<double> bernoulli;
    bernoulli.push_back(1);
    bernoulli.push_back(-0.5);
    double sum = 0;
    for (size_t i = 2; i <= n; ++i) {
        if (i % 2 == 1) {
            bernoulli.push_back(0);
            continue;
        } else {
            sum = 0;
            auto pascal = calc_pascal(i + 1);
            for (size_t k = 0; k < i; ++k) {
                double x = pascal[k + 2] * bernoulli[i - k - 1];
                sum += x;
            }
            bernoulli.push_back(-sum / (i + 1));
        }
    }
    return bernoulli;
}

int main() {
    std::vector<double> time;
    int n = 1000;
    for (int i = 1; i <= n; ++i) {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        auto bern = calc_bernoulli(i);
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        time.push_back(std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() / 1e9);
    }
    std::ofstream out("../data/bernoulli_time.txt");
    for (auto t : time) {
        out << t << std::endl;
    }
    out.close();
}