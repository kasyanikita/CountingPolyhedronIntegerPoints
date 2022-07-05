#include <iostream>
#include <string>
#include <vector>
#include <complex>
#include <algorithm>
#include <cmath>
#include <valarray>
#include "ToddPolynomial.h"

template<class TI, class TF>
class ToddFFT : public Todd<TI, TF> {
    typedef std::complex<TF> Complex;
    typedef std::valarray<Complex> CArray;
    void calc_todd() override;
    void ifft(CArray&);
    void fft(CArray &);
 public:
    ToddFFT(size_t, const std::vector<TI>&);
};

template<class TI, class TF>
ToddFFT<TI, TF>::ToddFFT(size_t m, const std::vector<TI>& xi): Todd<TI, TF>(m, xi) { }

template <class TI, class TF>
void ToddFFT<TI, TF>::fft(CArray &x) {
    size_t N = x.size(), k = N, n;
    TF thetaT = M_PI / N;
    Complex phiT = Complex(cos(thetaT), -sin(thetaT)), T;
    while (k > 1) {
        n = k;
        k >>= 1;
        phiT = phiT * phiT;
        T = 1.0;
        for (size_t l = 0; l < k; l++)
        {
            for (size_t a = l; a < N; a += n)
            {
                size_t b = a + k;
                Complex t = x[a] - x[b];
                x[a] += x[b];
                x[b] = t * T;
            }
            T *= phiT;
        }
    }
    size_t m = (size_t)log2(N);
    for (size_t a = 0; a < N; a++) {
        size_t b = a;
        b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
        b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
        b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
        b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
        b = ((b >> 16) | (b << 16)) >> (32 - m);
        if (b > a) {
            Complex t = x[a];
            x[a] = x[b];
            x[b] = t;
        }
    }
}

template <>
void ToddFFT<mpz_class, mpf_class>::fft(CArray &x) {
    size_t N = x.size(), k = N, n;
    mpf_class thetaT = M_PI / N;
    Complex phiT = Complex(cos(thetaT.get_d()), -sin(thetaT.get_d())), T;
    while (k > 1) {
        n = k;
        k >>= 1;
        phiT = phiT * phiT;
        T = 1.0;
        for (size_t l = 0; l < k; l++)
        {
            for (size_t a = l; a < N; a += n)
            {
                size_t b = a + k;
                Complex t = x[a] - x[b];
                x[a] += x[b];
                x[b] = t * T;
            }
            T *= phiT;
        }
    }
    size_t m = (size_t)log2(N);
    for (size_t a = 0; a < N; a++) {
        size_t b = a;
        b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
        b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
        b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
        b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
        b = ((b >> 16) | (b << 16)) >> (32 - m);
        if (b > a) {
            Complex t = x[a];
            x[a] = x[b];
            x[b] = t;
        }
    }
}


template <>
void ToddFFT<int64_t, mpf_class>::fft(CArray &x) {
    size_t N = x.size(), k = N, n;
    mpf_class thetaT = M_PI / N;
    Complex phiT = Complex(cos(thetaT.get_d()), -sin(thetaT.get_d())), T;
    while (k > 1) {
        n = k;
        k >>= 1;
        phiT = phiT * phiT;
        T = 1.0;
        for (size_t l = 0; l < k; l++)
        {
            for (size_t a = l; a < N; a += n)
            {
                size_t b = a + k;
                Complex t = x[a] - x[b];
                x[a] += x[b];
                x[b] = t * T;
            }
            T *= phiT;
        }
    }
    size_t m = (size_t)log2(N);
    for (size_t a = 0; a < N; a++) {
        size_t b = a;
        b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
        b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
        b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
        b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
        b = ((b >> 16) | (b << 16)) >> (32 - m);
        if (b > a) {
            Complex t = x[a];
            x[a] = x[b];
            x[b] = t;
        }
    }
}


template <class TI, class TF>
void ToddFFT<TI, TF>::ifft(CArray& x) {
    x = x.apply(std::conj);
 
    fft(x);
 
    x = x.apply(std::conj);
 
    // for (size_t i = 0; i < x.size(); ++i) {
    //     Complex den(x.size(), 0);
    //     x[i] = x[i] / den;
    // }
    x /= x.size();
}


template <class TI, class TF>
void ToddFFT<TI, TF>::calc_todd() {
    size_t n = 1;
    size_t deg = Todd<TI, TF>::todd.size() + Todd<TI, TF>::todd_part.size();
    while (n < deg) {
        n = n << 1;
    }
    CArray padX(n);
    CArray padY(n);
    for (size_t i = 0; i < Todd<TI, TF>::todd.size(); i++) {
        padX[i] = Todd<TI, TF>::todd[i];
		padY[i] = Todd<TI, TF>::todd_part[i];
    }
    fft(padX);
    fft(padY);
    padX *= padY;
    ifft(padX);
	for (size_t i = 0; i < Todd<TI, TF>::todd.size(); ++i) {
		Todd<TI, TF>::todd[i] = padX[i].real();
	}
}
