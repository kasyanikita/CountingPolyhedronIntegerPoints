#pragma once

#include <flint/fmpz_mat.h>
#include <vector>
#include <cassert>
#include <eigen3/Eigen/Dense>
#include <random>
#include <gmpxx.h>

namespace MSVP
{
    using namespace std;
    using namespace Eigen;

    using int_type = int64_t;
    using int_vec = Matrix<int_type, Dynamic, 1>;
    using int_mat = Matrix<int_type, Dynamic, Dynamic>;

    bool is_diagonal(const int_mat &A)
    {
        size_t n = A.cols();
        assert(n == A.rows());

        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                if (i != j && A(i, j) != 0)
                    return false;

        return true;
    }

    template <typename genT>
    int_mat generate_random_mat(size_t n, unsigned int interval, genT &gen)
    {
        uniform_int_distribution<int> distr(-interval, interval);

        int_mat res = int_mat::Zero(n, n);

        for (size_t i = 0; i < n; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                res(i, j) = distr(gen);
            }
        }

        return res;
    }

    namespace Internal
    {
        // void fill_flint_mat(const int_mat &A, fmpz_mat_t A_flint)
        // {
        //     size_t m = A.rows();
        //     size_t n = A.cols();

        //     fmpz_mat_init(A_flint, m, n);

        //     for (size_t i = 0; i < m; ++i)
        //         for (size_t j = 0; j < n; ++j)
        //             *fmpz_mat_entry(A_flint, i, j) = A(i, j);
        // }

        // int_mat from_flint_mat(size_t m, size_t n, fmpz_mat_t A_flint)
        // {
        //     int_mat res(m, n);

        //     for (size_t i = 0; i < m; ++i)
        //         for (size_t j = 0; j < n; ++j)
        //             res(i, j) = int_type(*fmpz_mat_entry(A_flint, i, j));

        //     return res;
        // }

        template <typename T>
        T euclid_bezout_nonnegative(
            const T &a,
            const T &b,
            T &u,
            T &v)
        {
            assert(a >= 0);
            assert(b >= 0);

            u = 1;
            v = 0;
            T u_new = v;
            T v_new = u;
            T r = a;
            T r_new = b;

            while (r_new != 0)
            {
                T q = r / r_new;

                T aux = u_new;
                u_new = u - q * u_new;
                u = aux;

                aux = v_new;
                v_new = v - q * v_new;
                v = aux;

                aux = r_new;
                r_new = r - q * r_new;
                r = aux;
            }

            return r;
        }

        template <typename T>
        T euclid_bezout(
            const T &a,
            const T &b,
            T &u,
            T &v)
        {
            T res = euclid_bezout_nonnegative<int_type>(abs(a), abs(b), u, v);

            if (a < 0)
                u = -u;
            if (b < 0)
                v = -v;

            return res;
        }

        template <typename vecT>
        size_t find_nonZero_index(const vecT &vec, size_t start)
        {
            size_t n = vec.size();

            for (size_t i = start; i < n; ++i)
            {
                if (vec[i] != 0)
                    return i;
            }

            return n;
        }

        bool make_swap(int_mat &A, int_mat &P, int_mat &Q, size_t pivot_i)
        {
            size_t n = A.cols();
            assert(n == A.rows());

            size_t non_zero_idx = find_nonZero_index(A.row(pivot_i), pivot_i + 1);
            if (non_zero_idx < n)
            {
                A.col(non_zero_idx).swap(A.col(pivot_i));
                Q.col(non_zero_idx).swap(Q.col(pivot_i));

                return true;
            }

            non_zero_idx = find_nonZero_index(A.col(pivot_i), pivot_i + 1);
            if (non_zero_idx < n)
            {
                A.row(non_zero_idx).swap(A.row(pivot_i));
                P.row(non_zero_idx).swap(P.row(pivot_i));

                return true;
            }

            return false;
        }

        template <typename vecT>
        size_t find_min_nonDivide(const vecT &vec, int_type d, size_t start)
        {
            size_t n = vec.size();
            int_type min = vec.array().abs().sum();
            size_t min_i = n;
            for (size_t i = start; i < n; ++i)
                if (vec[i] % d != 0)
                {
                    if (abs(vec[i]) < min)
                    {
                        min_i = i;
                        min = vec[i];
                    }
                }

            return min_i;
        }

        tuple<size_t, size_t> find_min_nonDivide(const int_mat &A, size_t pivot_i)
        {
            size_t n = A.cols();
            assert(n == A.rows());

            const int_type &d = A(pivot_i, pivot_i);
            assert(d != 0);

            size_t col_min_i = find_min_nonDivide(A.row(pivot_i), d, pivot_i + 1);
            size_t row_min_i = find_min_nonDivide(A.col(pivot_i), d, pivot_i + 1);

            if (col_min_i >= n)
                return make_tuple(row_min_i, n);

            if (row_min_i >= n)
                return make_tuple(n, col_min_i);

            if (abs(A(pivot_i, col_min_i)) < abs(A(row_min_i, pivot_i)))
                return make_tuple(n, col_min_i);

            return make_tuple(row_min_i, n);
        }
    }

    tuple<int_mat, int_mat, int_mat> my_DF(int_mat A)
    {
        size_t n = A.cols();
        assert(n == A.rows());

        int_mat P = int_mat::Identity(n, n),
                Q = int_mat::Identity(n, n);

        // pivot element
        size_t pivot_i = 0;
        while (pivot_i < n)
        {
            if (A(pivot_i, pivot_i) == 0)
            {
                // trying to find non-zero pivot
                // and make correpsonding swap
                bool succ = Internal::make_swap(A, P, Q, pivot_i);
                // if elelmts on the row and column with pivot_i index are all zeroes
                // then just increase pivot_i and move to the next iteration
                if (!succ)
                {
                    ++pivot_i;
                    continue;
                }
            }

            // if pivot is negative => make it positive
            if (A(pivot_i, pivot_i) < 0)
            {
                A.col(pivot_i) = -A.col(pivot_i);
                Q.col(pivot_i) = -Q.col(pivot_i);
            }

            // try to find a row or column index, such that
            // A(row_cand, pivot_i) or A(pivot_i, col_cand) is not multiple of A(pivot_i, pivot_i)
            // we try to find the minimal such element
            // if all elements are multiples of A(pivot_i, pivot_i)
            // then just row_cand == n, and col_cand == n
            size_t row_cand = 0, col_cand = 0;
            tie(row_cand, col_cand) = Internal::find_min_nonDivide(A, pivot_i);

            if (col_cand < n)
            {
                // make column extendet gcd computation
                int_type x = 0, y = 0;
                int_type gcd = Internal::euclid_bezout(A(pivot_i, pivot_i), A(pivot_i, col_cand), x, y);
                A.col(pivot_i) = (A.col(pivot_i) * x + A.col(col_cand) * y).eval();
                Q.col(pivot_i) = (Q.col(pivot_i) * x + Q.col(col_cand) * y).eval();
            }
            else if (row_cand < n)
            {
                // make row extendet gcd computation
                int_type x = 0, y = 0;
                Internal::euclid_bezout(A(pivot_i, pivot_i), A(row_cand, pivot_i), x, y);
                A.row(pivot_i) = (A.row(pivot_i) * x + A.row(row_cand) * y).eval();
                P.row(pivot_i) = (P.row(pivot_i) * x + P.row(row_cand) * y).eval();
            }
            else
            {
                // now, all elements are multiples of A(pivot_i, pivot_i)
                // lets make them 0
                // increase pivot_i
                // and move to the next iteration
                for (size_t col_i = pivot_i + 1; col_i < n; ++col_i)
                {
                    int_type div = A(pivot_i, col_i) / A(pivot_i, pivot_i);
                    assert(A(pivot_i, col_i) % A(pivot_i, pivot_i) == 0);

                    A.col(col_i) -= A.col(pivot_i) * div;
                    Q.col(col_i) -= Q.col(pivot_i) * div;
                }

                for (size_t row_i = pivot_i + 1; row_i < n; ++row_i)
                {
                    int_type div = A(row_i, pivot_i) / A(pivot_i, pivot_i);
                    assert(A(row_i, pivot_i) % A(pivot_i, pivot_i) == 0);

                    A(row_i, pivot_i) = 0;
                    P.row(row_i) -= P.row(pivot_i) * div;
                }

                ++pivot_i;
            }
        }

        return make_tuple(P, A, Q);
    }

    size_t unit_permut_normalize(int_mat &P, int_mat &S, int_mat &Q)
    {
        assert(is_diagonal(S));
        assert(P.rows() == P.cols() && S.rows() == S.cols() && Q.rows() == Q.cols());

        size_t n = S.cols();

        size_t nonZero_num = 0;
        size_t i = 0;
        while (i < n)
        {
            assert(S(i, i) > 0);
            if (S(i, i) > 1)
            {
                size_t swap_idx = n - nonZero_num - 1;
                swap(S(i, i), S(swap_idx, swap_idx));
                P.row(i).swap(P.row(swap_idx));
                Q.col(i).swap(Q.col(swap_idx));
                ++nonZero_num;
            }
            else
            {
                ++i;
            }
        }

        return nonZero_num;
    }
} // namespace MSVP
