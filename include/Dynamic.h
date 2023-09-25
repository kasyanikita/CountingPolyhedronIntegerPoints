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
        std::vector<std::vector<bool>> isComputed;
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
            for (int i = 0; i < n; ++i)
            {
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

        void normalize()
        {
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
            for (int i = 0; i < den.size(); ++i)
            {
                den[i] = den[i] / det;
            }

            // print_vector<int_t>(mat_vec_mult(Aadj, b),"Aadj * b");

            // normalize alpha
            for (int i = 0; i < exps.size(); ++i)
            {
                // std::cout << "Until exp: " << exps[i];
                exps[i] = (dot_product(c, mat_vec_mult(Aadj, b)) + exps[i]) / det;
                // std::cout << ", After exp: " << exps[i];
                // std::cout << std::endl;
            }

            d(n - 1, g[n]) = ExpPoly(exps, coeffs);
            // std::cout << d(n-1, g[n]) << std::endl;
        }

    public:
        Dynamic(std::vector<int_t> c_, std::vector<std::vector<int_t>> &A_, std::vector<int_t> &b_) : P(A_.size(), std::vector<int_t>(A_[0].size(), 0)),
                                                                                                      Q(A_.size(), std::vector<int_t>(A_[0].size(), 0)),
                                                                                                      S(A_.size(), 0)
        {
            c = c_;
            A = A_;
            b = b_;
            n = A.size();
        }

        void init()
        {
            get_SNF();
            g = calc_g(P, b, S);
            r = calc_r(g);
            h = calc_h(A);
            dp = init_dp();

            for (size_t i = 0; i < dp.size(); ++i)
            {
                isComputed.push_back(std::vector<bool>(dp[0].size(), false));
            }
            // print_matrix(P, "P");
            // print_vector<int_t>(S, "S");
            // print_vector<int_t>(g[0].get_components(), "g1");
            // print_vector<int_t>(g[1].get_components(), "g2");
            // print_vector<int_t>(g[2].get_components(), "g0");
            // print_vector<uint_t>(r, "r");
            // print_matrix(h, "h");
        }

        ExpPoly &d(uint_t k, GroupElement g)
        {
            return dp[k][g.get_idx()];
        }

        std::vector<std::vector<ExpPoly>> get_table()
        {
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

        std::vector<int_t> get_den()
        {
            return den;
        }

        ExpPoly get_final_poly()
        {
            return d(n - 1, g[n]);
        }

        void start()
        {
            init_first_layer();
            // std::cout <<  << std::endl;
            for (int k = 1; k < n; ++k)
            {
                for (int i = 0; i < dp[0].size(); ++i)
                {
                    // j = 0
                    // auto begin = std::chrono::high_resolution_clock::now();
                    auto ge = get_group_element_by_index(i, S);
                    // auto end = std::chrono::high_resolution_clock::now();
                    // group_element_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;
                    // std::cout << "Get group element by index: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;
                    // std::cout << std::endl;


                    // begin = std::chrono::high_resolution_clock::now();
                    auto scal_mul_vec = scalar_vector_mul(0, h[k]);
                    // end = std::chrono::high_resolution_clock::now();
                    // scalar_mult_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;
                    // std::cout << "Scalar mult: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;
                    // std::cout << std::endl;


                    // begin = std::chrono::high_resolution_clock::now();
                    auto dot_prod = -dot_product(c, scal_mul_vec);
                    // end = std::chrono::high_resolution_clock::now();
                    // dot_prod_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;
                    // std::cout << "Dot prod: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;
                    // std::cout << std::endl;


                    // begin = std::chrono::high_resolution_clock::now();
                    auto sum = d(k - 1, ge).monomial_multiply(dot_prod, 1);
                    // end = std::chrono::high_resolution_clock::now();
                    // monomial_mult_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;
                    // std::cout << "Monomial mult:" << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;
                    // std::cout << std::endl;


                    // j = 1 ... rk - 1
                    for (int j = 1; j < r[k]; ++j)
                    {
                        auto sum_part = d(k - 1, ge - j * g[k]).monomial_multiply(-dot_product(c, scalar_vector_mul(j, h[k])), 1);
                        // if (i == 0) {
                        //     std::cout << "f(k-1, g0 - " << j << "*" << "g2) = " << d(k - 1, ge - j * g[k]) << std::endl;
                        //     std::cout << "exp(-<c, >)"
                        // }


                        // begin = std::chrono::high_resolution_clock::now();
                        sum = sum + sum_part;
                        // end = std::chrono::high_resolution_clock::now();
                        // sum_part_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;
                        // std::cout << "Sum sum + sum part: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;
                        // std::cout << std::endl;
                    }

                    // save d(k, ge)
                    d(k, ge) = sum;
                }
            }
            // std::cout << d(n-1, g[n]) << std::endl;
            den = get_denominator();
            normalize();
            // std::cout << d(n - 1, g[n]) << std::endl;
        }

        void new_start()
        {
            auto res = (*this)(n - 1, g[n]);
            den = get_denominator();
            normalize();
        }

        ExpPoly operator()(uint_t k, const GroupElement &ge)
        {
            auto g_idx = ge.get_idx();
            if (isComputed[k][g_idx])
            {
                return dp[k][g_idx];
            }
            else
            {
                if (k != 0)
                {
                    auto sum = (*this)(k - 1, ge).monomial_multiply(-dot_product(c, scalar_vector_mul(0, h[k])), 1);

                    // j = 1 ... rk - 1
                    for (int j = 1; j < r[k]; ++j)
                    {
                        auto sum_part = (*this)(k - 1, ge - j * g[k]).monomial_multiply(-dot_product(c, scalar_vector_mul(j, h[k])), 1);
                        // if (i == 0) {
                        //     std::cout << "f(k-1, g0 - " << j << "*" << "g2) = " << d(k - 1, ge - j * g[k]) << std::endl;
                        //     std::cout << "exp(-<c, >)"
                        // }
                        sum = sum + sum_part;
                    }

                    // save d(k, ge)
                    dp[k][g_idx] = sum;
                    isComputed[k][g_idx] = true;
                }
                else
                {
                    // k = 0
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
                    d(k, ge).init({exp}, {coeff});

                    isComputed[k][g_idx] = true;
                }
            }
            return dp[k][g_idx];
        }
    };
}