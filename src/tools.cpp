#include "Tools.h"

namespace GroupIP {

std::vector<GroupIP::GroupElement> CalculateGroupElements(const Matrix& P, const Vector &snf_diagonal,
                                                          const Vector &b) {

  int n = b.size();
  std::vector<GroupElement> g(n + 1, GroupElement(snf_diagonal));
  g[n].assign(mat_vec_mult(P, b));

  for (int i = 0; i < n; ++i) {
    g[i].assign(get_mat_col(P, i));
  }

  return g;
}

Matrix CalculateH(const Matrix &A) {
  Matrix h(A.size(), Vector(A.size()));
  auto Aadj = CalculateAdjugateMatrix(A);

  for (int i = 0; i < A.size(); ++i) {
    for (int j = 0; j < A[i].size(); ++j) {
      h[j][i] = Aadj[i][j];
    }
  }
  
  return h;
}

ExpPoly vertex_normalize(const ExpPoly &num, std::vector<int_t> &den,
                      const Matrix &A, const Vector &b, const Vector &c) {
  auto num_poly = num.get_poly();

  // init exps and coeffs
  std::vector<ExpPoly::coeff_t> coeffs;
  std::vector<ExpPoly::exp_t> exps;
  for (const auto &[e, c] : num_poly) {
    coeffs.push_back(c);
    exps.push_back(e);
  }

  // get adjugate matrix and det of the matrix A
  auto Aadj = CalculateAdjugateMatrix(A);
  auto det = get_determinant(A);

  // normalize denominator
  for (int i = 0; i < den.size(); ++i) {
    den[i] = den[i] / det;
  }

  // normalize alpha
  for (int i = 0; i < exps.size(); ++i) {
    exps[i] = (dot_product(c, mat_vec_mult(Aadj, b)) + exps[i]) / det;
  }

  return ExpPoly(exps, coeffs);
}

std::vector<int_t> get_denominator(const Vector &c,
                                   std::vector<GroupElement> &g,
                                   const std::vector<std::vector<int_t>> &h) {
  std::vector<int_t> res;
  for (int i = 0; i < h.size(); ++i) {
    int_t sum = 0;
    for (int j = 0; j < h[0].size(); ++j) {
      sum += c[j] * h[i][j];
    }
    res.push_back((g[i].getOrder() * sum));
  }
  return res;
}

template <class T>
void print_vector(std::vector<T> v, std::string name) {
  std::cout << name << ": ";
  for (int i = 0; i < v.size(); ++i) {
    std::cout << v[i] << " ";
  }
  std::cout << '\n';
}

void print_matrix(std::vector<std::vector<int_t>> &M, std::string name) {
  std::cout << name << ":\n";
  for (int i = 0; i < M.size(); ++i) {
    for (int j = 0; j < M[i].size(); ++j) {
      std::cout << M[i][j] << " ";
    }
    std::cout << '\n';
  }
  std::cout << "\n";
}

std::vector<int_t> scalar_vector_mul(int_t s, const std::vector<int_t> &v) {
  std::vector<int_t> res(v.size());
  for (int i = 0; i < v.size(); ++i) {
    res[i] = s * v[i];
  }
  return res;
}

int_t dot_product(const std::vector<int_t> &a, const std::vector<int_t> &b) {
  int_t res = 0;
  for (int i = 0; i < a.size(); ++i) {
    res += a[i] * b[i];
  }
  return res;
}

std::vector<int_t> get_mat_col(const std::vector<std::vector<int_t>> &M,
                               int k) {
  std::vector<int_t> res;
  for (int i = 0; i < M.size(); ++i) {
    res.push_back(M[i][k]);
  }
  return res;
}

std::vector<int_t> mat_vec_mult(const std::vector<std::vector<int_t>> &M,
                                const std::vector<int_t> &v) {
  std::vector<int_t> res;
  for (int i = 0; i < M.size(); ++i) {
    int_t sum = 0;
    for (int j = 0; j < M[i].size(); ++j) {
      sum += M[i][j] * v[j];
    }
    res.push_back(sum);
  }
  return res;
}

int_t calc_s(const GroupElement &g, const GroupElement &g0, uint_t r0) {
  int_t res = -1;
  for (int i = 0; i < r0; ++i) {
    if (i * g0 == g) {
      res = i;
      break;
    }
  }
  return res;
}

GroupElement get_group_element_by_index(int_t idx, std::vector<int_t> &S) {
  int_t div = 1;
  for (int i = 0; i < S.size() - 1; ++i) {
    div *= S[i];
  }
  std::vector<int_t> comp(S.size(), 0);
  for (int i = 0; i < S.size() - 1; ++i) {
    comp[S.size() - i - 1] = idx / div;
    idx = idx % div;
    div /= S[S.size() - i - 2];
  }
  comp[0] = idx;
  GroupElement res(S);
  res.assign(comp);
  return res;
}

}  // namespace GroupIP

int_t get_random_number(int_t min, int_t max) {
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_int_distribution<int_t> dist(min, max);
  return dist(rng);
}

int_t get_determinant(const GroupIP::Matrix &A) {
  fmpz_mat_t Af;
  fmpz_t det;
  fmpz_mat_init(Af, A.size(), A[0].size());
  for (int i = 0; i < A.size(); ++i) {
    for (int j = 0; j < A[i].size(); ++j) {
      auto val = fmpz_mat_entry(Af, i, j);
      *val = A[i][j];
    }
  }
  fmpz_mat_det(det, Af);
  return *det;
}

std::vector<std::vector<int_t>> get_sub_matrix(
    std::vector<std::vector<int_t>> &A, int_t row_to_exclude) {
  std::vector<std::vector<int_t>> res;
  for (int i = 0; i < A.size(); ++i) {
    if (i != row_to_exclude) res.push_back(A[i]);
  }
  return res;
}

std::vector<int_t> get_sub_vector(std::vector<int_t> &b, int_t row_to_exclude) {
  std::vector<int_t> res;
  for (int i = 0; i < b.size(); ++i) {
    if (i != row_to_exclude) res.push_back(b[i]);
  }
  return res;
}

std::vector<int_t> gen_rand_vector(int_t n, int_t a, int_t b) {
  std::vector<int_t> v(n);
  for (int j = 0; j < n; ++j) {
    v[j] = get_random_number(a, b);
  }

  return v;
}

Matrix CalculateAdjugateMatrix(const Matrix &A) {
  fmpz_mat_t Aadj;
  fmpz_t den;
  fmpz_t det;
  std::vector<std::vector<int_t>> res(A.size(),
                                      std::vector<int_t>(A[0].size()));
  fmpz_mat_init(Aadj, A.size(), A[0].size());
  for (int i = 0; i < A.size(); ++i) {
    for (int j = 0; j < A[i].size(); ++j) {
      auto val = fmpz_mat_entry(Aadj, i, j);
      *val = A[i][j];
    }
  }

  fmpz_mat_det(det, Aadj);
  int is_not_singular = fmpz_mat_inv(Aadj, den, Aadj);
  if (is_not_singular == 0) {
    throw std::domain_error("Matrix A is singular");
  }

  if (*det != *den) {
    std::cout << "det != den" << std::endl;
  }
  fmpz_divexact(den, det, den);
  fmpz_mat_scalar_divexact_fmpz(Aadj, Aadj, den);

  for (int i = 0; i < A.size(); ++i) {
    for (int j = 0; j < A[i].size(); ++j) {
      auto val = fmpz_mat_entry(Aadj, i, j);
      res[i][j] = fmpz_get_d(val);
    }
  }
  return res;
}

std::vector<int_t> adj_c_multiply(std::vector<std::vector<int_t>> &Aadj,
                                  std::vector<int_t> &c) {
  std::vector<int_t> res(Aadj[0].size(), 0);
  for (int i = 0; i < Aadj[0].size(); ++i) {
    for (int j = 0; j < Aadj.size(); ++j) {
      res[i] = res[i] + c[j] * Aadj[j][i];
    }
  }
  return res;
}

bool is_any_zero(std::vector<int_t> &v) {
  for (int i = 0; i < v.size(); ++i) {
    if (v[i] == 0) return true;
  }
  return false;
}

bool check_sub(std::vector<std::vector<int_t>> &A_sub, std::vector<int_t> &c) {
  auto A_sub_adj = CalculateAdjugateMatrix(A_sub);
  auto res = adj_c_multiply(A_sub_adj, c);
  bool flag = is_any_zero(res);
  return flag;
}

bool check_subs(std::vector<std::vector<std::vector<int_t>>> &A_subs,
                std::vector<int_t> &c) {
  for (int i = 0; i < A_subs.size(); ++i) {
    if (check_sub(A_subs[i], c)) {
      return true;
    }
  }
  return false;
}

std::vector<mpf_class> get_sum_exps(
    int_t n, std::vector<ExpPoly::exp_t> exps,  // std::vector<mpf_class>
    std::vector<ExpPoly::coeff_t> coeffs, std::vector<mpf_class> todd) {
  std::vector<std::vector<int_t>> powers(exps.size(),
                                         std::vector<int_t>(n + 1, 1));
  fill_powers(powers, exps);
  // print_matrix(powers, "Powers");
  auto fact = get_factorial(n);

  std::vector<mpf_class> summa;
  for (int i = 0; i < exps.size(); ++i) {
    mpf_class sum = 0;
    for (int j = 0; j <= n; ++j) {
      sum += (powers[i][j] * todd[n - j]) / fact[j];
    }
    // std::cout << "exp: " << exps[i] << ", sum: " << sum << std::endl;
    summa.push_back(coeffs[i] * sum);
  }
  return summa;
}

std::vector<int_t> get_factorial(int_t n) {
    std::vector<int_t> fact(n + 1, 1);
    for (int i = 2; i <= n; ++i) {
        fact[i] = fact[i - 1] * i;
    }
    return fact;
}

void fill_powers(std::vector<std::vector<int_t>>& powers,
std::vector<ExpPoly::exp_t> exps) {
    for (int i = 0; i < exps.size(); ++i) {
        for (int j = 1; j < powers[i].size(); ++j) {
            powers[i][j] = powers[i][j - 1] * exps[i];
        }
    }
}