#include <iostream>
#include <vector>
#include <random>
#include <flint/fmpz_mat.h>
#include "SmithNormalForm.h"

std::vector<std::vector<int>> gen_rand_mat(int, int, int);
bool test(int, int, int);
void run_tests(int, int, int, int);
bool check_PAQ(fmpz_mat_t, std::vector<std::vector<int>>&, std::vector<std::vector<int>>&,
                 std::vector<std::vector<int>>&);
void convert_vec_to_fmpz_mat(fmpz_mat_t, const std::vector<std::vector<int>> &);

int main()
{
    int tries = 1000;
    int n = 5;
    int a = 0;
    int b = 5;
    run_tests(tries, n, a, b);
}

void run_tests(int tries, int n, int a, int b) {
    int success = 0;
    for (int i = 0; i < tries; ++i) {
        if (test(n, a, b)) {
            ++success;
        } else {
            std::cout << "Error!\n";
        }
    }
    std::cout << success << "/" << tries << '\n';
}

bool test(int n, int a, int b) {
    auto A = gen_rand_mat(n, a, b);
    std::vector<std::vector<int>> A_copy(n, std::vector<int>(n));

    // init matrices
    std::vector<std::vector<int>> S(n, std::vector<int>(n));
    std::vector<std::vector<int>> P(n, std::vector<int>(n, 0));
    std::vector<std::vector<int>> Q(n, std::vector<int>(n, 0));
    fmpz_mat_t Sf;
    fmpz_mat_t Af;
    fmpz_mat_init(Af, n, n);
    fmpz_mat_init(Sf, n, n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            auto val = fmpz_mat_entry(Af, i, j);
            *val = A[i][j];
            A_copy[i][j] = A[i][j];
        }
    }

    // calculate snf
    fmpz_mat_snf(Sf, Af);
    SNF(A, P, Q);

    // convert fmpz_mat to std::vector
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            auto val = fmpz_mat_entry(Sf, i, j);
            S[i][j] = *val;
        }
    }

    bool res = S == A;
    res = res && check_PAQ(Sf, P, A_copy, Q);

    return res;
}

bool check_PAQ(fmpz_mat_t S, std::vector<std::vector<int>> & P, std::vector<std::vector<int>> & A,
             std::vector<std::vector<int>> & Q)
{
    fmpz_mat_t Pf;
    fmpz_mat_t Af;
    fmpz_mat_t Qf;
    fmpz_mat_init(Pf, P.size(), P[0].size());
    fmpz_mat_init(Af, A.size(), A[0].size());
    fmpz_mat_init(Qf, Q.size(), Q[0].size());

    convert_vec_to_fmpz_mat(Pf, P);
    convert_vec_to_fmpz_mat(Af, A);
    convert_vec_to_fmpz_mat(Qf, Q);

    fmpz_mat_mul(Pf, Pf, Af);
    fmpz_mat_mul(Pf, Pf, Qf);

    if (fmpz_mat_equal(S, Pf) > 0) return true;
    return false;
}

std::vector<std::vector<int>> gen_rand_mat(int n, int a, int b)
{
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist(a, b);

    std::vector<std::vector<int>> A(n, std::vector<int>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A[i][j] = dist(rng);
        }
    }
    
    return A;
}

void convert_vec_to_fmpz_mat(fmpz_mat_t Mf, const std::vector<std::vector<int>> & M) {
    for (int i = 0; i < M.size(); ++i) {
        for (int j = 0; j < M[i].size(); ++j) {
            auto val = fmpz_mat_entry(Mf, i, j);
            *val = M[i][j];
        }
    }
}