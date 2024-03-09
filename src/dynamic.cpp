#include "Dynamic.h"

namespace GroupIP {

std::vector<std::vector<ExpPoly>> Dynamic::init_dp() {
  auto S = snf_.GetDiagonalOfS();
  int_t n_cols = 1;
  for (int i = 0; i < n; ++i) {
    n_cols *= S[i];
  }
  std::vector<std::vector<ExpPoly>> res(n, std::vector<ExpPoly>(n_cols));
  return res;
}

Dynamic::Dynamic(const std::vector<std::vector<int_t>> &A,
                 const std::vector<int_t> &b, const std::vector<int_t> &c)
    : snf_(A), b_(b), c_(c), aux_(A, b, snf_) {}

void Dynamic::init() {
  snf_.CalculateSNF();
  n = snf_.GetDiagonalOfS().size();
  aux_.Init();
  dp = init_dp();

  for (size_t i = 0; i < dp.size(); ++i) {
    isComputed.push_back(std::vector<bool>(dp[0].size(), false));
  }
}

ExpPoly &Dynamic::d(uint_t k, GroupElement g) { return dp[k][g.get_idx()]; }

ExpPoly Dynamic::operator()(uint_t k, const GroupElement &ge) {
  auto h = aux_.GetH();
  auto g = aux_.GetG();
  auto g_idx = ge.get_idx();

  if (isComputed[k][g_idx]) {
    return dp[k][g_idx];
  } else {
    if (k != 0) {
      auto sum = (*this)(k - 1, ge).monomial_multiply(
          -dot_product(c_, scalar_vector_mul(0, h[k])), 1);

      // j = 1 ... rk - 1
      for (int j = 1; j < g[k].getOrder(); ++j) {
        auto sum_part = (*this)(k - 1, ge - j * g[k])
                            .monomial_multiply(
                                -dot_product(c_, scalar_vector_mul(j, h[k])), 1);
        sum = sum + sum_part;
      }

      // save d(k, ge)
      dp[k][g_idx] = sum;
      isComputed[k][g_idx] = true;
    } else {
      // k = 0
      auto s = calc_s(ge, g[0], g[0].getOrder());
      int_t exp = 0;
      uint_t coeff = 0;
      if (s != -1) {
        exp = -dot_product(c_, scalar_vector_mul(s, h[0]));
        coeff = 1;
      }
      d(k, ge).init({exp}, {coeff});
      isComputed[k][g_idx] = true;
    }
  }
  return dp[k][g_idx];
}
}  // namespace GroupIP