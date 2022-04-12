#include "ToddPolynomial.h"

std::vector<uint64_t> Todd::calc_pascal(int n) {
    // Get Pascal's triangle n-th row

    std::vector<uint64_t> res(n + 1, 1);
    int i;
    for (i = 0; i < n; ++i) {
        res[i + 1] = res[i] * (n - i) / (i + 1);
    }
    return res;
}

void Todd::calc_bernulli(int n) {
    // Calculate bernulli numbers

    bernulli.push_back(1);
    bernulli.push_back(-0.5);
    int i;
    int k;
    double sum;
    for (i = 2; i <= n; ++i) {
        if (i % 2 == 1) {
            bernulli.push_back(0);
            continue;
        }
        sum = 0;
        auto pascal = calc_pascal(i + 1);
        for (k = 0; k < i; ++k) {
            sum += pascal[k + 2] * bernulli[i - k - 1];
        }
        bernulli.push_back((-1.0 / (i + 1)) * sum);
    }
}

Todd::Todd(uint64_t m, std::vector<double> xi) {
    calc_bernulli(m);
}

double Todd::get_bernulli(int i) {
    return bernulli[i];
}
