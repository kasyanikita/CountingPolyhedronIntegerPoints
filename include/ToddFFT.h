#ifndef INCLUDE_TODDFFT_H_
#define INCLUDE_TODDFFT_H_

#include <complex>
#include <algorithm>
#include <cmath>
#include <valarray>
#include <flint/fmpz_poly.h>
#include <flint/fmpq_poly.h>
#include <gmp.h>
#include <gmpxx.h>
#include "ToddPolynomial.h"

template <class TI, class TF>
class ToddFFT : public Todd<TI, TF>
{
    std::vector<mpf_class> res;
    void calc_todd() override;

public:
    ToddFFT(size_t, const std::vector<TI> &);
    std::vector<mpf_class> &get_todd_mpf();
};

template <class TI, class TF>
ToddFFT<TI, TF>::ToddFFT(size_t m, const std::vector<TI> &xi) : Todd<TI, TF>(m, xi) {}

template <class TI, class TF>
void ToddFFT<TI, TF>::calc_todd()
{
    fmpq_poly_t a, b, res;
    size_t size = Todd<TI, TF>::todd.size();
    fmpq_poly_init2(a, size);
    fmpq_poly_init2(b, size);
    fmpq_poly_init2(res, size);
    for (size_t i = 0; i < size; ++i)
    {
        fmpq_poly_set_coeff_mpq(a, i, mpq_class(Todd<TI, TF>::todd[i]).get_mpq_t());
        fmpq_poly_set_coeff_mpq(b, i, mpq_class(Todd<TI, TF>::todd_part[i]).get_mpq_t());
    }
    fmpq_poly_mullow(res, a, b, size);
    mpq_t tmp;
    mpq_init(tmp);
    for (size_t i = 0; i < size; ++i)
    {
        fmpq_poly_get_coeff_mpq(tmp, res, i);
        Todd<TI, TF>::todd[i] = mpq_class(tmp);
    }
}

template <class TI, class TF>
std::vector<mpf_class> &ToddFFT<TI, TF>::get_todd_mpf()
{
    for (size_t i = 0; i < Todd<TI, TF>::todd.size(); ++i)
    {
        res.push_back(mpf_class(Todd<TI, TF>::todd[i]));
    }
    return res;
}

#endif // INCLUDE_TODDFFT_H_
