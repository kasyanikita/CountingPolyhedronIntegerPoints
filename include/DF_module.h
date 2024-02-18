#pragma once

// #include <vector>
#include <cassert>
#include <eigen3/Eigen/Dense>
// #include <random>
// #include <gmpxx.h>
#include <iostream>

using namespace Eigen;

namespace Internal
{
    template <typename T>
    T euclid_bezout_nonnegative(
        const T &a,
        const T &b,
        T &u, T &v,
        T &u0, T &v0)
    {
        assert(a >= 0);
        assert(b >= 0);

        u = 1;
        u0 = 0;
        v = 0;
        v0 = 1;

        T r = a;
        T r0 = b;

        while (r0 != 0)
        {
            T q = r / r0;

            T aux = u0;
            u0 = u - q * u0;
            u = aux;

            aux = v0;
            v0 = v - q * v0;
            v = aux;

            aux = r0;
            r0 = r - q * r0;
            r = aux;
        }

        return r;
    }

    template <typename T>
    T euclid_bezout(
        const T &a,
        const T &b,
        T &u,
        T &v,
        T &u0,
        T &v0)
    {
        T res = euclid_bezout_nonnegative<T>(abs(a), abs(b), u, v, u0, v0);

        if (a < 0)
        {
            u = -u;
            u0 = -u0;
        }
        if (b < 0)
        {
            v = -v;
            v0 = -v0;
        }

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
}

template <typename intT>
class my_diagonal_form_algo
{
public:
    using int_type = intT;
    using int_vec = Eigen::Matrix<int_type, Eigen::Dynamic, 1>;
    using int_mat = Eigen::Matrix<int_type, Eigen::Dynamic, Eigen::Dynamic>;

private:
    bool make_swap(int_mat &A, int_mat &P, int_mat &Q, size_t pivot_i)
    {
        size_t n = A.cols();
        assert(n == A.rows());

        size_t non_zero_idx = Internal::find_nonZero_index(A.row(pivot_i), pivot_i + 1);
        if (non_zero_idx < n)
        {
            A.col(non_zero_idx).swap(A.col(pivot_i));
            Q.col(non_zero_idx).swap(Q.col(pivot_i));

            return true;
        }

        non_zero_idx = Internal::find_nonZero_index(A.col(pivot_i), pivot_i + 1);
        if (non_zero_idx < n)
        {
            A.row(non_zero_idx).swap(A.row(pivot_i));
            P.row(non_zero_idx).swap(P.row(pivot_i));

            return true;
        }

        return false;
    }

    size_t find_min_nonDivide(const int_vec &vec, int_type d, size_t start)
    {
        size_t n = vec.size();
        int_type min = vec.array().abs().maxCoeff() + 1;
        size_t min_i = n;
        for (size_t i = start; i < n; ++i)
            if (vec[i] % d != 0)
            {
                if (abs(vec[i]) < min)
                {
                    min_i = i;
                    min = abs(vec[i]);
                }
            }

        return min_i;
    }

    std::tuple<size_t, size_t> find_min_nonDivide(const int_mat &A, size_t pivot_i)
    {
        size_t n = A.cols();
        assert(n == A.rows());

        const int_type &d = A(pivot_i, pivot_i);
        assert(d != 0);

        size_t col_min_i = find_min_nonDivide(A.row(pivot_i), d, pivot_i + 1);
        size_t row_min_i = find_min_nonDivide(A.col(pivot_i), d, pivot_i + 1);

        if (col_min_i >= n)
            return std::make_tuple(row_min_i, n);

        if (row_min_i >= n)
            return std::make_tuple(n, col_min_i);

        if (abs(A(pivot_i, col_min_i)) < abs(A(row_min_i, pivot_i)))
            return std::make_tuple(n, col_min_i);

        return std::make_tuple(row_min_i, n);
    }

public:
    std::tuple<int_mat, int_mat, int_mat> compute(int_mat A)
    {
        // std::cout << A << endl;
        size_t n = A.cols();
        assert(n == A.rows());

        int_mat P = int_mat::Identity(n, n),
                Q = int_mat::Identity(n, n);

        // pivot element
        size_t pivot_i = 0;
        while (pivot_i + 1 < n)
        {
            if (A(pivot_i, pivot_i) == 0)
            {
                // trying to find non-zero pivot
                // and make correpsonding swap
                bool succ = make_swap(A, P, Q, pivot_i);
                // if elelmts on the row and column with pivot_i index are all zeroes
                // then just increase pivot_i and move to the next iteration
                if (!succ)
                {
                    ++pivot_i;
                    continue;
                }
            }
            // std::cout << A << endl;

            // if pivot is negative => make it positive
            if (A(pivot_i, pivot_i) < 0)
            {
                A.col(pivot_i) = -A.col(pivot_i);
                Q.col(pivot_i) = -Q.col(pivot_i);
            }
            // std::cout << A << endl;

            // try to find a row or column index, such that
            // A(row_cand, pivot_i) or A(pivot_i, col_cand) is not multiple of A(pivot_i, pivot_i)
            // we try to find the minimal such element
            // if all elements are multiples of A(pivot_i, pivot_i)
            // then just row_cand == n, and col_cand == n
            size_t row_cand = 0, col_cand = 0;
            std::tie(row_cand, col_cand) = find_min_nonDivide(A, pivot_i);

            if (col_cand < n)
            {
                // make column extendet gcd computation
                int_type x1 = 0, x2 = 0, y1 = 0, y2 = 0;
                int_type gcd = Internal::euclid_bezout(A(pivot_i, pivot_i), A(pivot_i, col_cand),
                                                       x1, y1, x2, y2);

                int_vec col1 = A.col(pivot_i) * x1 + A.col(col_cand) * y1;
                int_vec col2 = A.col(pivot_i) * x2 + A.col(col_cand) * y2;
                A.col(pivot_i) = col1;
                A.col(col_cand) = col2;

                col1 = Q.col(pivot_i) * x1 + Q.col(col_cand) * y1;
                col2 = Q.col(pivot_i) * x2 + Q.col(col_cand) * y2;
                Q.col(pivot_i) = col1;
                Q.col(col_cand) = col2;

                assert(A(pivot_i, pivot_i) == gcd);
                assert(A(pivot_i, col_cand) == 0);

                // std::cout << A << endl;
            }
            else if (row_cand < n)
            {
                // make row extendet gcd computation
                int_type x1 = 0, x2 = 0, y1 = 0, y2 = 0;
                int_type gcd = Internal::euclid_bezout(A(pivot_i, pivot_i), A(row_cand, pivot_i),
                                                       x1, y1, x2, y2);
                int_vec row1 = A.row(pivot_i) * x1 + A.row(row_cand) * y1;
                int_vec row2 = A.row(pivot_i) * x2 + A.row(row_cand) * y2;

                A.row(pivot_i) = row1;
                A.row(row_cand) = row2;

                row1 = P.row(pivot_i) * x1 + P.row(row_cand) * y1;
                row2 = P.row(pivot_i) * x2 + P.row(row_cand) * y2;
                P.row(pivot_i) = row1;
                P.row(row_cand) = row2;

                assert(A(pivot_i, pivot_i) == gcd);
                assert(A(row_cand, pivot_i) == 0);

                // std::cout << A << endl;
            }
            else
            {
                // std::cout << A << endl;
                //  now, all elements are multiples of A(pivot_i, pivot_i)
                //  lets make them 0
                //  increase pivot_i
                //  and move to the next iteration
                for (size_t col_i = pivot_i + 1; col_i < n; ++col_i)
                {
                    int_type div = A(pivot_i, col_i) / A(pivot_i, pivot_i);
                    // if (A(pivot_i, col_i) % A(pivot_i, pivot_i) != 0)
                    // {
                    //     std::cout << A(pivot_i, pivot_i) << ", " << A(pivot_i, col_i) << std::endl;
                    //     std::cout << A << std::endl;
                    //     std::cin.ignore();
                    // }
                    assert(A(pivot_i, col_i) % A(pivot_i, pivot_i) == 0);

                    A.col(col_i) -= A.col(pivot_i) * div;
                    Q.col(col_i) -= Q.col(pivot_i) * div;
                }

                for (size_t row_i = pivot_i + 1; row_i < n; ++row_i)
                {
                    int_type div = A(row_i, pivot_i) / A(pivot_i, pivot_i);
                    // if (A(row_i, pivot_i) % A(pivot_i, pivot_i) != 0)
                    // {
                    //     std::cout << A(pivot_i, pivot_i) << ", " << A(row_i, pivot_i) << std::endl;
                    //     std::cout << A << std::endl;
                    //     std::cin.ignore();
                    // }
                    assert(A(row_i, pivot_i) % A(pivot_i, pivot_i) == 0);

                    A(row_i, pivot_i) = 0;
                    P.row(row_i) -= P.row(pivot_i) * div;
                }

                ++pivot_i;
            }
        }

        return std::make_tuple(P, A, Q);
    }

    size_t normalize(int_mat &P, int_mat &S, int_mat &Q)
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
};
