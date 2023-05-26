#ifndef _TOOLS_H_
#define _TOOLS_H_

#include <stdexcept>
#include <vector>
#include <random>
#include <flint/fmpz_mat.h>
#include "global_defs.h"

using namespace GroupIP;

int_t get_random_number(int_t min, int_t max)
{
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<int_t> dist(min, max);
    return dist(rng);
}

int_t get_determinant(std::vector<std::vector<int_t>> &A) {
    fmpz_mat_t Af;
    fmpz_t det;
    fmpz_mat_init(Af, A.size(), A[0].size());
    for (int i = 0; i < A.size(); ++i)
    {
        for (int j = 0; j < A[i].size(); ++j)
        {
            auto val = fmpz_mat_entry(Af, i, j);
            *val = A[i][j];
        }
    }
    fmpz_mat_det(det, Af);
    return *det;
}

std::vector<std::vector<int_t>> get_sub_matrix(std::vector<std::vector<int_t>> &A, int_t row_to_exclude) {
    std::vector<std::vector<int_t>> res;
    for (int i = 0; i < A.size(); ++i) {
        if (i != row_to_exclude) res.push_back(A[i]);
    }
    return res;
}

std::vector<int_t> get_sub_vector(std::vector<int_t> &b, int_t row_to_exclude)
{
    std::vector<int_t> res;
    for (int i = 0; i < b.size(); ++i)
    {
        if (i != row_to_exclude)
            res.push_back(b[i]);
    }
    return res;
}

std::vector<int_t> gen_rand_vector(int_t n, int_t a, int_t b) {
    std::vector<int_t> v(n);
    for (int j = 0; j < n; ++j)
    {
        v[j] = get_random_number(a, b);
    }

    return v;
}

std::vector<std::vector<int_t>> get_adjugate(const std::vector<std::vector<int_t>>& A) {
    fmpz_mat_t Aadj;
    fmpz_t den;
    fmpz_t det;
    std::vector<std::vector<int_t>> res(A.size(), std::vector<int_t>(A[0].size()));
    fmpz_mat_init(Aadj, A.size(), A[0].size());
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A[i].size(); ++j) {
            auto val = fmpz_mat_entry(Aadj, i, j);
            *val = A[i][j];
        }
    }

    fmpz_mat_det(det, Aadj);
    int is_not_singular = fmpz_mat_inv(Aadj, den, Aadj);
    if (is_not_singular == 0)
    {
        throw std::domain_error("Matrix A is singular");
    }

    if (*det != *den)
    {
        std::cout << "det != den" << std::endl;
    }
    fmpz_divexact(den, det, den);
    fmpz_mat_scalar_divexact_fmpz(Aadj, Aadj, den);

    for (int i = 0; i < A.size(); ++i)
    {
        for (int j = 0; j < A[i].size(); ++j)
        {
            auto val = fmpz_mat_entry(Aadj, i, j);
            res[i][j] = fmpz_get_d(val);
        }
    }
    return res;
}

std::vector<int_t> adj_c_multiply(std::vector<std::vector<int_t>>& Aadj, std::vector<int_t>& c) {
    std::vector<int_t> res(Aadj[0].size(), 0);
    for (int i = 0; i < Aadj[0].size(); ++i)
    {
        for (int j = 0; j < Aadj.size(); ++j)
        {
            res[i] = res[i] + c[j] * Aadj[j][i];
        }
    }
    return res;
}

bool is_any_zero(std::vector<int_t>& v) {
    for (int i = 0; i < v.size(); ++i) {
        if (v[i] == 0) return true;
    }
    return false;
}

bool check_sub(std::vector<std::vector<int_t>> &A_sub, std::vector<int_t> &c) {
    auto A_sub_adj = get_adjugate(A_sub);
    // print_matrix(A_sub_adj, "A_adj");
    auto res = adj_c_multiply(A_sub_adj, c);
    // print_vector<int_t>(res, "c * Aadj");
    return (!is_any_zero(res));
}

bool check_subs(std::vector<std::vector<std::vector<int_t>>> &A_subs, std::vector<int_t> &c)
{
    for (int i = 0; i < A_subs.size(); ++i) {
        if (!check_sub(A_subs[i], c)) {
            return false;
        }
    }
    return true;
}

std::vector<int_t> get_c_vector(std::vector<std::vector<int_t>> &A, int_t alpha)
{
    bool flag = false;
    std::vector<int_t> c;
    std::vector<std::vector<std::vector<int_t>>> A_subs;
    for (int i = 0; i < A.size(); ++i) {
        auto A_sub = get_sub_matrix(A, i);
        A_subs.push_back(A_sub);
    }
    while (!flag){
        c = gen_rand_vector(A_subs[0].size(), -alpha, alpha);
        flag = check_subs(A_subs, c);
        // print_vector<int_t>(c, "c");
    }
    return c;
}

#endif