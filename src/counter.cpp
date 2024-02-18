#include "counter.h"

Counter::Counter(Matrix& A, Vector& b) {
  A_ = A;
  b_ = b; 
}

mpf_class Counter::count_simple_cone(Matrix &A, Vector &b, Vector &c) {

  //init dynamic
  Dynamic d(A, b, c);
  int_t n = A.size();

  d.init();
  std::cout << "Init success!\n";
  d.new_start();
  std::cout << "New start success!\n";
  auto res = d.get_final_poly();
  auto den = d.get_den();
  std::vector<mpz_class> den_mpz(den.size());

  for (int i = 0; i < den.size(); ++i) {
    den_mpz[i] = mpz_class(den[i]);
  }

  // get coeffs and exps
  auto poly = res.get_poly();
  std::vector<ExpPoly::coeff_t> coeffs;
  std::vector<ExpPoly::exp_t> exps;
  for (const auto &[exp, coeff] : poly) {
    exps.push_back(exp);
    coeffs.push_back(coeff);
  }

  ToddPoly_FFT<mpz_class, mpf_class> todd_poly(n, den_mpz);
  todd_poly.init();
  auto todd = todd_poly.get_todd();
  auto sum_alpha = get_sum_exps(n, exps, coeffs, todd);
  mpf_class numer = 0;
  mpz_class denom = 1;

  for (int i = 0; i < sum_alpha.size(); ++i) {
    numer += sum_alpha[i];
  }

  for (int i = 0; i < den_mpz.size(); ++i) {
    denom *= den_mpz[i];
  }

  return numer / denom;
}

mpf_class Counter::count_integer_points(Matrix& A, Vector& b) {
  mpf_class res = 0;
  // auto det = get_determinant(A);
  // auto c = get_c_vector(A, A.size());
  std::vector<int_t> c = {1, -1};
  for (int i = 0; i < A.size(); ++i) {
    auto Asub = get_sub_matrix(A, i);
    auto bsub = get_sub_vector(b, i);
    res += count_simple_cone(Asub, bsub, c);
  }
  return res;
}