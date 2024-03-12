#include "counter.h"

mpf_class PolyEvaluation(Matrix &A, Vector &b, Vector &c) {
  int_t n = A.size();
  SNFClass snf(A);
  snf.CalculateSNF();

  auto snf_diagonal = snf.GetDiagonalOfS();
  auto P = snf.GetP();
  auto g = CalculateGroupElements(P, snf_diagonal, b);
  auto h = CalculateH(A);

  Dynamic d(c, g, h);
  d.Init(snf_diagonal);

  auto num = d(n - 1, g[n]);
  auto den = GroupIP::get_denominator(c, g, h);

  num = GroupIP::vertex_normalize(num, den, A, b, c);

  std::vector<mpz_class> den_mpz(den.size());

  for (int i = 0; i < den.size(); ++i) {
    den_mpz[i] = mpz_class(den[i]);
  }

  // get coeffs and exps
  auto poly = num.get_poly();
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
