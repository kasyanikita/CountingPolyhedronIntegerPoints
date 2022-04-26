#ifndef INCLUDE_TODDPOLYNOMIAL_H_
#define INCLUDE_TODDPOLYNOMIAL_H_

#include <vector>
#include <gmpxx.h>

template <class TI, class TF>
class Todd {
 private:
    std::vector<TF> xi;
    std::vector<TF> bernulli;
    std::vector<TF> todd;
    std::vector<TI> calc_pascal(uint64_t) const;
    void calc_bernulli(uint64_t);
    void calc_todd(uint64_t, const std::vector<TF>&);
 public:
    Todd(uint64_t, const std::vector<TF>&);
    TF get_bernulli(uint64_t) const;
    TF get_todd(uint64_t) const;
};


template <class TI, class TF>
std::vector<TI> Todd<TI, TF>::calc_pascal(uint64_t n) const {
    // Get Pascal's triangle n-th row

    std::vector<TI> res(n + 1, 1);
    uint64_t i;
    for (i = 0; i < n; ++i) {
        res[i + 1] = res[i] * (n - i) / (i + 1);
    }
    return res;
}

template <class TI, class TF>
void Todd<TI, TF>::calc_bernulli(uint64_t n) {
    // Calculate bernulli numbers

    bernulli.push_back(1);
    bernulli.push_back(-0.5);
    uint64_t i;
    uint64_t k;
    TF sum;
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

template <class TI, class TF>
void Todd<TI, TF>::calc_todd(uint64_t m, const std::vector<TF>& xi) {
    std::vector<TF> next_todd_part(m + 1);
    std::vector<TF> next_todd(m + 1);
    TI fact = 1;
    TF curr_pow_xi = 1;
    TF next_pow_xi = 1;
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
            TF sum = 0;
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

template <class TI, class TF>
Todd<TI, TF>::Todd(uint64_t m, const std::vector<TF>& xi): todd(m + 1) {
   calc_bernulli(m);
   calc_todd(m, xi);
}

template <class TI, class TF>
TF Todd<TI, TF>::get_bernulli(uint64_t i) const {
   return bernulli[i];
}

template <class TI, class TF>
TF Todd<TI, TF>::get_todd(uint64_t i) const {
   return todd[i];
}

#endif  // INCLUDE_TODDPOLYNOMIAL_H_
