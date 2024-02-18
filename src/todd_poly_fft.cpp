#include "ToddPoly_FFT.h"

template <class TI, class TF>
ToddPoly_FFT<TI, TF>::ToddPoly_FFT(GroupIP::int_t m, const std::vector<TI> &xi)
    : ToddPoly<TI, TF>(m, xi) {}

template <class TI, class TF>
void ToddPoly_FFT<TI, TF>::calc_todd() {
  fmpq_poly_t a, b, res;
  size_t size = ToddPoly<TI, TF>::todd.size();
  fmpq_poly_init2(a, size);
  fmpq_poly_init2(b, size);
  fmpq_poly_init2(res, size);
  for (size_t i = 0; i < size; ++i) {
    fmpq_poly_set_coeff_mpq(a, i,
                            mpq_class(ToddPoly<TI, TF>::todd[i]).get_mpq_t());
    fmpq_poly_set_coeff_mpq(
        b, i, mpq_class(ToddPoly<TI, TF>::todd_part[i]).get_mpq_t());
  }
  fmpq_poly_mullow(res, a, b, size);
  mpq_t tmp;
  mpq_init(tmp);
  for (size_t i = 0; i < size; ++i) {
    fmpq_poly_get_coeff_mpq(tmp, res, i);
    ToddPoly<TI, TF>::todd[i] = mpq_class(tmp);
  }
}

template <class TI, class TF>
std::vector<mpf_class> &ToddPoly_FFT<TI, TF>::get_todd_mpf() {
  for (size_t i = 0; i < ToddPoly<TI, TF>::todd.size(); ++i) {
    res.push_back(mpf_class(ToddPoly<TI, TF>::todd[i]));
  }
  return res;
}