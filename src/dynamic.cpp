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

std::vector<int_t> Dynamic::get_denominator() {
  auto h = aux_.GetH();
  auto g = aux_.GetG();
  std::vector<int_t> res;
  for (int i = 0; i < h.size(); ++i) {
    int_t sum = 0;
    for (int j = 0; j < h[0].size(); ++j) {
      sum += c_[j] * h[i][j];
    }
    res.push_back((g[i].getOrder() * sum));
  }
  return res;
}

void Dynamic::normalize() {
  auto A = snf_.GetA();
  auto g = aux_.GetG();
  auto res = d(n - 1, g[n]);
  auto res_poly = res.get_poly();

  // init exps and coeffs
  std::vector<ExpPoly::coeff_t> coeffs;
  std::vector<ExpPoly::exp_t> exps;
  for (const auto &[e, c] : res_poly) {
    coeffs.push_back(c);
    exps.push_back(e);
  }

  // get adjugate materix and det of the matrix A
  auto Aadj = get_adjugate(A);
  auto det = get_determinant(A);

  // normalize denominator
  for (int i = 0; i < den.size(); ++i) {
    den[i] = den[i] / det;
  }

  // normalize alpha
  for (int i = 0; i < exps.size(); ++i) {
    exps[i] = (dot_product(c_, mat_vec_mult(Aadj, b_)) + exps[i]) / det;
  }

  d(n - 1, g[n]) = ExpPoly(exps, coeffs);
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

std::vector<std::vector<ExpPoly>> Dynamic::get_table() { return dp; }

void Dynamic::init_first_layer() {
  auto S = snf_.GetDiagonalOfS();
  auto g = aux_.GetG();
  auto h = aux_.GetH();
  
  for (int i = 0; i < dp[0].size(); ++i) {
    auto ge = get_group_element_by_index(i, S);
    auto s = calc_s(ge, g[0], g[0].getOrder());

    int_t exp = 0;
    uint_t coeff = 0;

    if (s != -1) {
      exp = -dot_product(c_, scalar_vector_mul(s, h[0]));
      coeff = 1;
    }
    d(0, ge).init({exp}, {coeff});
  }
}

std::vector<int_t> Dynamic::get_den() { return den; }

ExpPoly Dynamic::get_final_poly() {
  auto g = aux_.GetG();
  return d(n - 1, g[n]);
}

// void Dynamic::start() {
//   init_first_layer();
//   // std::cout <<  << std::endl;
//   for (int k = 1; k < n; ++k) {
//     for (int i = 0; i < dp[0].size(); ++i) {

//       auto ge = get_group_element_by_index(i, S);
//       auto scal_mul_vec = scalar_vector_mul(0, h[k]);
//       auto dot_prod = -dot_product(c, scal_mul_vec);
//       auto sum = d(k - 1, ge).monomial_multiply(dot_prod, 1);

//       // j = 1 ... rk - 1
//       for (int j = 1; j < r[k]; ++j) {
//         auto sum_part = d(k - 1, ge - j * g[k])
//                             .monomial_multiply(
//                                 -dot_product(c, scalar_vector_mul(j, h[k])),
//                                 1);

//         sum = sum + sum_part;
//       }

//       // save d(k, ge)
//       d(k, ge) = sum;
//     }
//   }
//   den = get_denominator();
//   normalize();
// }

void Dynamic::new_start() {
  auto g = aux_.GetG();
  auto res = (*this)(n - 1, g[n]);
  den = get_denominator();
  normalize();
}

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