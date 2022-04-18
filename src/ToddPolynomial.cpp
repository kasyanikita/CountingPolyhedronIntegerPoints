#include "ToddPolynomial.h"

std::vector<uint64_t> Todd::calc_pascal(int n) const {
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

void Todd::calc_todd(uint64_t m, const std::vector<double>& xi) {
    std::vector<double> next_todd_part(m + 1);
    std::vector<double> next_todd(m + 1);
    uint64_t fact = 1;
    double curr_pow_xi = 1;
    double next_pow_xi = 1;
    uint64_t i;
    uint64_t j;
    uint64_t k;
    for (i = 0; i <= m; ++i) {
        todd[i] = curr_pow_xi * bernulli[i] / fact;
        if (xi.size() != 1) {
            next_todd_part[i] = next_pow_xi * bernulli[i] / fact;
            next_pow_xi *= -xi[1];
        }
        if (i > 0) fact *= (i + 1);
        curr_pow_xi *= -xi[0];
    }

    for (i = 1; i < xi.size(); ++i) {
        fact = 1;
        next_pow_xi = 1;
        for (j = 0; j <= m; ++j) {
            double sum = 0;
            for (k = 0; k <= j; ++k) {
                sum += todd[j - k] * next_todd_part[k];
            }
            next_todd[j] = sum;
        }

        if (i != xi.size() - 1) {
            for (j = 0; j <= m; ++j) {
                next_todd_part[j] = next_pow_xi * bernulli[j] / fact;
                next_pow_xi *= -xi[i + 1];
                if (j > 0) fact *= (j + 1);
            }
        }
        todd.swap(next_todd);
    }
}

Todd::Todd(uint64_t m, const std::vector<double>& xi): todd(m + 1) {
    calc_bernulli(m);
    calc_todd(m, xi);
}

double Todd::get_bernulli(int i) const {
    return bernulli[i];
}

double Todd::get_todd(uint64_t i) const {
    return todd[i];
}
