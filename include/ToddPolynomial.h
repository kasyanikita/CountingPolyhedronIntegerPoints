#ifndef INCLUDE_TODDPOLYNOMIAL_H_
#define INCLUDE_TODDPOLYNOMIAL_H_

#include <vector>
#include <iostream>
#include <gmpxx.h>

template <class TI, class TF>
class Todd {
 protected:
    size_t m;
    std::vector<TI> xi;
    std::vector<TF> bernoulli;
    std::vector<TF> todd;
    std::vector<TF> todd_part;
    std::vector<TF> tmp_todd;
    std::vector<TI> calc_pascal(size_t) const;
    void calc_bernoulli();
    void update_part(TI);
    void init_todd();
    void virtual calc_todd();
    void calc_todd_arr();
 public:
    Todd(size_t, const std::vector<TI>&);
    void init();
    const std::vector<TF>& get_bernoulli() const;
    const std::vector<TF>& get_todd() const;
    std::vector<TF>& get_bernoulli();
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
void Todd<TI, TF>::calc_bernoulli() {
    // Calculate bernoulli numbers from 0 to m
    bernoulli.push_back(1);
    bernoulli.push_back(-0.5);
    TF sum = 0;
    for (size_t i = 2; i <= m; ++i) {
        if (i % 2 == 1) {
            bernoulli.push_back(0);
            continue;
        } else {
            sum = 0;
            auto pascal = calc_pascal(i + 1);
            for (size_t k = 0; k < i; ++k) {
                TF x = pascal[k + 2] * bernoulli[i - k - 1];
                sum += x;
            }
            bernoulli.push_back(-sum / (i + 1));
        }
    }
}


template <>
void Todd<mpz_class, mpf_class>::calc_bernoulli() {
    // Calculate bernoulli numbers from 0 to m
    bernoulli.push_back(1);
    bernoulli.push_back(-0.5);
    mpf_class sum(0, 500);
    for (size_t i = 2; i <= m; ++i) {
        if (i % 2 == 1) {
            bernoulli.push_back(0);
            continue;
        } else {
            sum = 0;
            auto pascal = calc_pascal(i + 1);
            for (size_t k = 0; k < i; ++k) {
                sum += pascal[k + 2] * bernoulli[i - k - 1];
            }
            bernoulli.push_back(-sum / (i + 1));
        }
    }
}


template<class TI, class TF>
void Todd<TI, TF>::update_part(TI x) {
    todd_part[0] = 1;
    TI fact = 1;
    TI pow_x = 1;
    for (size_t i = 1; i <= m; ++i) {
        fact *= i;
        pow_x *= -x;
        todd_part[i] = pow_x * bernoulli[i] / fact;
    }
}


template<class TI, class TF>
void Todd<TI, TF>::init_todd() {
    todd = todd_part;
}


template<class TI, class TF>
void Todd<TI, TF>::calc_todd() {
    for (size_t i = 0; i <= m; ++i) {
        tmp_todd[i] = 0;
        for (size_t j = 0; j <= i; ++j) {
            tmp_todd[i] += todd[i - j] * todd_part[j];
        }
    }
    todd.swap(tmp_todd);
}


template <class TI, class TF>
void Todd<TI, TF>::calc_todd_arr() {
    // Calculate todd polynomials of degree from 0 to m
    update_part(xi[0]);
    init_todd();
    for (size_t i = 1; i < xi.size(); ++i) {
        update_part(xi[i]);
        calc_todd();
    }
}


template <class TI, class TF>
Todd<TI, TF>::Todd(size_t _m, const std::vector<TI>& _xi): todd(_m + 1), todd_part(_m + 1), tmp_todd(_m + 1) {
    m = _m;
    xi = _xi;
}


template <class TI, class TF>
void Todd<TI, TF>::init() {
    calc_bernoulli();
    calc_todd_arr();
}


template <class TI, class TF>
const std::vector<TF>& Todd<TI, TF>::get_bernoulli() const {
    return bernoulli;
}


template <class TI, class TF>
std::vector<TF>& Todd<TI, TF>::get_bernoulli() {
    return bernoulli;
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
