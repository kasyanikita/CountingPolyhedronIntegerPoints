#pragma once

#include <stdexcept>
#include "global_defs.h"
#include "GroupElement.h"
#include "ExpPoly.h"
#include "SmithNormalForm.h"
#include "Tools.h"
#include <flint/fmpz_mat.h>

namespace GroupIP
{

    std::vector<int_t> scalar_vector_mul(int_t s, const std::vector<int_t>& v)
    {
        std::vector<int_t> res(v.size());
        // std::cout << "s: " << s << " vector: ";
        for (int i = 0; i < v.size(); ++i)
        {
            // std::cout << v[i] << " ";
            res[i] = s * v[i];
        }
        // std::cout << std::endl;
        return res;
    }

    int_t dot_product(const std::vector<int_t>& a, const std::vector<int_t>& b)
    {
        int_t res = 0;
        // std::cout << "Dot product: ";
        for (int i = 0; i < a.size(); ++i)
        {
            // std::cout << a[i] << "*" << b[i] << " ";
            res += a[i] * b[i];
            // if (i != a.size() - 1) std::cout << "+ ";
        }
        // std::cout << std::endl;
        return res;
    }

    std::vector<int_t> get_mat_col(const std::vector<std::vector<int_t>> &M, int k)
    {
        std::vector<int_t> res;
        for (int i = 0; i < M.size(); ++i)
        {
            res.push_back(M[i][k]);
        }
        return res;
    }

    std::vector<int_t> mat_vec_mult(const std::vector<std::vector<int_t>> &M, const std::vector<int_t> &v)
    {
        std::vector<int_t> res;
        for (int i = 0; i < M.size(); ++i)
        {
            int_t sum = 0;
            for (int j = 0; j < M[i].size(); ++j)
            {
                sum += M[i][j] * v[j];
            }
            res.push_back(sum);
        }
        return res;
    }

    std::vector<GroupElement> calc_g(const std::vector<std::vector<int_t>> &P,
                                     const std::vector<int_t> &b, const std::vector<int_t> &S)
    {
        int n = b.size();
        std::vector<GroupElement> g_vec(n + 1, GroupElement(S));
        g_vec[n].assign(mat_vec_mult(P, b));
        for (int i = 0; i < n; ++i)
        {
            g_vec[i].assign(get_mat_col(P, i));
        }
        return g_vec;
    }

    std::vector<uint_t> calc_r(const std::vector<GroupElement> &g_vec)
    {
        std::vector<int_t> e(g_vec.size() - 1, 0);
        std::vector<uint_t> r_vec;
        for (int i = 0; i < g_vec.size(); ++i)
        {
            int r = 1;
            auto sum = g_vec[i];
            while (sum.get_components() != e)
            {
                sum += g_vec[i];
                ++r;
            }
            r_vec.push_back(r);
        }
        return r_vec;
    }

    int_t calc_s(GroupElement& g, GroupElement& g0, uint_t r0)
    {
        int_t res = -1;
        for (int i = 0; i < r0; ++i) {
            if (i * g0 == g) {
                res = i;
                break;
            } 
        }
        return res;
    }

    std::vector<std::vector<int_t>> calc_h(const std::vector<std::vector<int_t>> &A)
    {
        // init variables
        std::vector<std::vector<int_t>> h(A.size(), std::vector<int_t>(A.size()));
        fmpz_mat_t A_flint;
        fmpz_t den;
        fmpz_t det;
        fmpz_mat_t Ainv;
        fmpz_mat_init(Ainv, A.size(), A[0].size());
        fmpz_mat_init(A_flint, A.size(), A[0].size());

        // std::vector to fmpz_mat_t
        for (int i = 0; i < A.size(); ++i)
        {
            for (int j = 0; j < A[i].size(); ++j)
            {
                auto val = fmpz_mat_entry(A_flint, i, j);
                *val = A[i][j];
            }
        }

        // get inverse matrix
        fmpz_mat_det(det, A_flint);
        int res = fmpz_mat_inv(Ainv, den, A_flint);
        if (res == 0)
        {
            throw std::domain_error("Matrix A is singular");
        }

        // calculate adjugate matrix
        if (*det != *den)
        {
            std::cout << "det != den" << std::endl;
        }
        fmpz_divexact(den, det, den);
        fmpz_mat_scalar_divexact_fmpz(Ainv, Ainv, den);

        // fmpz_mat_t to std::vector
        for (int i = 0; i < A.size(); ++i)
        {
            for (int j = 0; j < A[i].size(); ++j)
            {
                auto val = fmpz_mat_entry(Ainv, i, j);
                h[j][i] = fmpz_get_d(val);
            }
        }
        return h;
    }

    GroupElement get_group_element_by_index(int_t idx, std::vector<int_t>& S) {
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

    class Dynamic
    {
        int_t n;
        std::vector<int_t> c;
        std::vector<std::vector<int_t>> A;
        std::vector<std::vector<int_t>> P;
        std::vector<std::vector<int_t>> Q;
        std::vector<int_t> S;
        std::vector<int_t> b;
        std::vector<GroupElement> g;
        std::vector<std::vector<int_t>> h;
        std::vector<uint_t> r;
        std::vector<std::vector<ExpPoly>> dp;
        std::vector<int_t> den;

        void get_SNF()
        {
            std::vector<std::vector<int_t>> S_mat = A;
            SNF(S_mat, P, Q);
            for (int i = 0; i < n; ++i)
            {
                S[i] = S_mat[i][i];
            }
        }

        std::vector<std::vector<ExpPoly>> init_dp()
        {
            int_t n_cols = 1;
            for (int i = 0; i < n; ++i) {
                n_cols *= S[i];
            }
            std::vector<std::vector<ExpPoly>> res(n, std::vector<ExpPoly>(n_cols));
            return res;
        }

        std::vector<int_t> get_denominator()
        {
            std::vector<int_t> res;
            for (int i = 0; i < h.size(); ++i)
            {
                int_t sum = 0;
                for (int j = 0; j < h[0].size(); ++j)
                {
                    sum += c[j] * h[i][j];
                }
                res.push_back((r[i] * sum));
            }
            return res;
        }

        void normalize() {
            auto res = d(n - 1, g[n]);
            auto res_poly = res.get_poly();

            // init exps and coeffs
            std::vector<ExpPoly::coeff_t> coeffs;
            std::vector<ExpPoly::exp_t> exps;
            for (const auto &[e, c] : res_poly)
            {
                coeffs.push_back(c);
                exps.push_back(e);
            }

            // get adjugate materix and det of the matrix A
            auto Aadj = get_adjugate(A);
            auto det = get_determinant(A);
            // std::cout << "Determinant:" << det << std::endl;
            // print_matrix(Aadj, "Adjugate");
            // print_vector<int_t>(b, "b");

            // normalize denominator
            for (int i = 0; i < den.size(); ++i) {
                den[i] = den[i] / det;
            }

            // print_vector<int_t>(mat_vec_mult(Aadj, b),"Aadj * b");

            // normalize alpha
            for (int i = 0; i < exps.size(); ++i) {
                // std::cout << "Until exp: " << exps[i];
                exps[i] = (dot_product(c, mat_vec_mult(Aadj, b)) + exps[i]) / det;
                // std::cout << ", After exp: " << exps[i];
                // std::cout << std::endl;
            }

            d(n-1, g[n]) = ExpPoly(exps, coeffs);
            // std::cout << d(n-1, g[n]) << std::endl;
        }

    public:
        Dynamic(std::vector<int_t> c_, std::vector<std::vector<int_t>> &A_, std::vector<int_t> &b_) :
         P(A_.size(), std::vector<int_t>(A_[0].size(), 0)),
         Q(A_.size(), std::vector<int_t>(A_[0].size(), 0)),
         S(A_.size(), 0)
        {
            c = c_;
            A = A_;
            b = b_;
            n = A.size();
        }

        void init() {
            get_SNF();
            g = calc_g(P, b, S);
            r = calc_r(g);
            h = calc_h(A);
            dp = init_dp();
            // print_matrix(P, "P");
            // print_vector<int_t>(S, "S");
            // print_vector<int_t>(g[0].get_components(), "g1");
            // print_vector<int_t>(g[1].get_components(), "g2");
            // print_vector<int_t>(g[2].get_components(), "g0");
            // print_vector<uint_t>(r, "r");
            // print_matrix(h, "h");
        }

        ExpPoly& d(uint_t k, GroupElement g) {
            return dp[k][g.get_idx()];
        }

        std::vector<std::vector<ExpPoly>> get_table() {
            return dp;
        }

        void init_first_layer()
        {
            // std::cout << dp[0].size() << std::endl;
            for (int i = 0; i < dp[0].size(); ++i)
            {
                auto ge = get_group_element_by_index(i, S);
                auto s = calc_s(ge, g[0], r[0]);
                // std::cout << "min(s*g1=gi) = " << s << std::endl;
                // std::cout << s << std::endl;
                int_t exp = 0;
                uint_t coeff = 0;
                // std::cout << "Init vector: ";
                // for (int i = 0; i < ge.get_components().size(); ++i) {
                //     std::cout << ge.get_components()[i] << " ";
                // }
                // std::cout << '\n';
                if (s != -1)
                {
                    exp = -dot_product(c, scalar_vector_mul(s, h[0]));
                    coeff = 1;
                }
                d(0, ge).init({exp}, {coeff});
                // std::cout << d(0, ge) << std::endl;
                // std::cout << "=================\n\n";
            }
        }

        std::vector<int_t> get_den() {
            return den;
        }

        ExpPoly get_final_poly() { 
            return d(n-1, g[n]);
        }

        void start() {
            init_first_layer();
            // std::cout <<  << std::endl;
            for (int k = 1; k < n; ++k) {
                for (int i = 0; i < dp[0].size(); ++i) {
                    // j = 0
                    auto ge = get_group_element_by_index(i, S);
                    auto sum = d(k - 1, ge).monomial_multiply(-dot_product(c, scalar_vector_mul(0, h[k])), 1);

                    // j = 1 ... rk - 1
                    for (int j = 1; j < r[k]; ++j) {
                        auto sum_part = d(k - 1, ge - j * g[k]).monomial_multiply(-dot_product(c, scalar_vector_mul(j, h[k])), 1);
                        // if (i == 0) {
                        //     std::cout << "f(k-1, g0 - " << j << "*" << "g2) = " << d(k - 1, ge - j * g[k]) << std::endl;
                        //     std::cout << "exp(-<c, >)"
                        // }
                        sum = sum + sum_part;
                    }

                    // save d(k, ge)
                    d(k, ge) = sum;
                }
            }
            den = get_denominator();
            normalize();
            // std::cout << d(n-1, g[n]) << std::endl;
        }
    };
}