#ifndef INCLUDE_TODDPOLYNOMIAL_H_
#define INCLUDE_TODDPOLYNOMIAL_H_

#include <vector>
#include <iostream>
#include <gmpxx.h>

template <class TI, class TF>
class Todd {
 private:
    size_t m;
    std::vector<TI> xi;
    std::vector<TF> bernulli;
    std::vector<TF> todd;
    std::vector<TI> calc_pascal(size_t) const;
    void calc_bernulli();
    void calc_todd();
 public:
    Todd(size_t, const std::vector<TI>&);
    void init();
    const std::vector<TF>& get_bernulli() const;
    const std::vector<TF>& get_todd() const;
    std::vector<TF>& get_bernulli();
    std::vector<TF>& get_todd();
};


template <class TI, class TF>
std::vector<TI> Todd<TI, TF>::calc_pascal(size_t n) const {
    // Get Pascal's triangle n-th row
    std::vector<TI> res(n + 1, 1);
    for (size_t i = 0; i < n; ++i) {
        res[i + 1] = res[i] * (n - i) / (i + 1);
    }
    return res;
}

template <class TI, class TF>
void Todd<TI, TF>::calc_bernulli() {
    // Calculate bernulli numbers from 0 to m
    bernulli.push_back(1);
    bernulli.push_back(-0.5);
    TF sum = 0;
    for (size_t i = 2; i <= m; ++i) {
        if (i % 2 == 1) {
            bernulli.push_back(0);
            continue;
        } else {
            sum = 0;
            auto pascal = calc_pascal(i + 1);
            for (size_t k = 0; k < i; ++k) {
                sum += pascal[k + 2] * bernulli[i - k - 1];
            }
            bernulli.push_back(-sum / (i + 1));
        }
    }
}

template <>
void Todd<mpz_class, mpf_class>::calc_bernulli() {
    // Calculate bernulli numbers from 0 to m
    bernulli.push_back(1);
    bernulli.push_back(-0.5);
    mpf_class sum(0, 128);
    for (size_t i = 2; i <= m; ++i) {
        if (i % 2 == 1) {
            bernulli.push_back(0);
            continue;
        } else {
            sum = 0;
            auto pascal = calc_pascal(i + 1);
            for (size_t k = 0; k < i; ++k) {
                sum += pascal[k + 2] * bernulli[i - k - 1];
            }
            bernulli.push_back(-sum / (i + 1));
        }
    }
}

template <class TI, class TF>
void Todd<TI, TF>::calc_todd() {
    // Calculate todd polynomial of degree from 0 to m
    std::vector<TF> next_todd_part(m + 1);
    std::vector<TF> next_todd(m + 1);
    TI fact = 1;
    TI curr_pow_xi = 1;
    TI next_pow_xi = 1;

    for (size_t i = 0; i <= m; ++i) {
        todd[i] = curr_pow_xi * bernulli[i] / fact;
        if (xi.size() != 1) {
            next_todd_part[i] = next_pow_xi * bernulli[i] / fact;
            next_pow_xi *= -xi[1];
        }
        if (i > 0) fact *= (i + 1);
        curr_pow_xi *= -xi[0];
    }

    for (size_t i = 1; i < xi.size(); ++i) {
        fact = 1;
        next_pow_xi = 1;
        for (size_t j = 0; j <= m; ++j) {
            TF sum = 0;
            for (size_t k = 0; k <= j; ++k) {
                sum += todd[j - k] * next_todd_part[k];
            }
            next_todd[j] = sum;
        }

        if (i != xi.size() - 1) {
            for (size_t j = 0; j <= m; ++j) {
                next_todd_part[j] = next_pow_xi * bernulli[j] / fact;
                next_pow_xi *= -xi[i + 1];
                if (j > 0) fact *= (j + 1);
            }
        }
        todd.swap(next_todd);
    }
}

template <class TI, class TF>
Todd<TI, TF>::Todd(size_t _m, const std::vector<TI>& _xi): todd(_m + 1) {
    m = _m;
    xi = _xi;
}

template <class TI, class TF>
void Todd<TI, TF>::init() {
    calc_bernulli();
    calc_todd();
}

template <class TI, class TF>
const std::vector<TF>& Todd<TI, TF>::get_bernulli() const {
    return bernulli;
}

template <class TI, class TF>
std::vector<TF>& Todd<TI, TF>::get_bernulli() {
    return bernulli;
}

template <class TI, class TF>
const std::vector<TF>& Todd<TI, TF>::get_todd() const {
    return todd;
}

template <class TI, class TF>
std::vector<TF>& Todd<TI, TF>::get_todd() {
    return todd;
}

#endif  // INCLUDE_TODDPOLYNOMIAL_H_
