#include "ToddPoly.h"

template <>
void ToddPoly<mpz_class, mpf_class>::calc_bernoulli() {
  // Calculate bernoulli numbers from 0 to m
  bernoulli.push_back(1);
  bernoulli.push_back(-0.5);
  mpf_class sum(0, 500);
  for (size_t i = 2; i <= m; ++i) {
    if (i % 2 == 1) {
      bernoulli.push_back(0);
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